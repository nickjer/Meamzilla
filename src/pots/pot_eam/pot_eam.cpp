#include "pot_eam.h"
#include "atom_eam.h"
#include "pair_eam.h"
#include "../../input.h"
#include "../../config.h"
#include "../../cell.h"
#include "../../atom_vec.h"
#include "../../error.h"
#include "../../pot_list.h"
#include "../../style_basis.h"

#include <fstream>
#include <iomanip>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

PotEAM::PotEAM(class Meamzilla *mmz, int ntypes) : PotPair(mmz, ntypes), embed_extrap_(0)
{
  // Double count pairs
  dbl_cnt_pair_ = 1;

  // Resize list of potentials
  npot_fns_ = 3;
  pot_fns_.resize(npot_fns_);

  // Density potential (rho) has n potentials (depending on atom j, the neighbor) and is a radial potential
  // For types A and B, you have rho_A, rho_B
  pot_fns_[1].setup_pot_fns(ntypes, PotFns::AlloyType::Atom_j, PotFns::FnType::Radial);

  // Embedding fn has n potentials (depending on atom i, the original atom) and is NOT a radial potential
  // For types A and B, you have F_A, F_B
  pot_fns_[2].setup_pot_fns(ntypes, PotFns::AlloyType::Atom_i, PotFns::FnType::Undefined);
}

/* ---------------------------------------------------------------------- */

PotEAM::~PotEAM()
{
  // dtor
}

/* ----------------------------------------------------------------------
   clone potential
------------------------------------------------------------------------- */

PotEAM* PotEAM::clone() const
{
  return new PotEAM(*this);
}

/* ----------------------------------------------------------------------
   initialize potential
------------------------------------------------------------------------- */

void PotEAM::init()
{
  if (mmz->input) {
    mmz->input->parse("embed_extrap",    0, embed_extrap_);
  }

  return;
}

/* ----------------------------------------------------------------------
   check that a potential is similar enough to this
   potential for pairs/triplets to have identical info
------------------------------------------------------------------------- */

void PotEAM::check_if_similar(Potential *pot_ptr) const
{
  PotEAM& pot = *static_cast<PotEAM *>(pot_ptr);

  // Compare phi
  pot_fns_[0].check_if_similar(pot.pot_fns_[0]);  // error will be thrown if not similar

  // Compare rho
  pot_fns_[1].check_if_similar(pot.pot_fns_[1]);  // error will be thrown if not similar

  // Do not compare embedding fn

  return;
}

/* ----------------------------------------------------------------------
   trap non-root procs in infinite loop
------------------------------------------------------------------------- */

void PotEAM::compute_trap(const Comm& comm, int flag)
{
  Vector<double> tmp_vec; // nothing is ever written to this variable

  if (comm.is_root() && flag)
    is_trapped_ = 1;
  else if (comm.is_root())
    is_trapped_ = 0;

  do {
    comm.bcast(&flag, 1, MPI_INT, comm.get_root());

    if (flag == 2)
      compute(comm);
    else if (flag == 3)
      compute(comm, &tmp_vec);
    else if (flag == 4)
      compute_densities(comm, &tmp_vec);

  } while (flag && !comm.is_root());

  return;
}

/* ----------------------------------------------------------------------
   compute error sum/vector
------------------------------------------------------------------------- */

double PotEAM::compute(const Comm& comm, ErrorVec *error_vec)
{
  // Untrap procs if necessary
  if (is_trapped_) {
    if (error_vec) {
      int flag = 3;
      comm.bcast(&flag, 1, MPI_INT, comm.get_root());
    } else {
      int flag = 2;
      comm.bcast(&flag, 1, MPI_INT, comm.get_root());
    }
  }

  // Initialize potential on all procs
  initialize_pot(comm, error_vec);

  // Initialize potential by resetting forces
  initialize_compute(comm);

  // Setup local variables for potential functions
  std::vector<Basis *> phi_fns;
  for (Basis*& fn : pot_fns_[0].fns)
    phi_fns.push_back(fn);
  std::vector<Basis *> rho_fns;
  for (Basis*& fn : pot_fns_[1].fns)
    rho_fns.push_back(fn);
  std::vector<Basis *> F_fns;
  for (Basis*& fn : pot_fns_[2].fns)
    F_fns.push_back(fn);

  // Set up constraint error (error from density going out of bounds of embedding function)
  Vector<double> constraint_err(mmz->config->ncells,0.0);

  // Loop over all atoms in atomvec
  for (Atom*& atom_i_ptr : mmz->atomvec->atoms) {
    // Make temporary atom and cell
    AtomEAM &atom_i = *(static_cast<AtomEAM *>(atom_i_ptr));
    Cell &cell = *mmz->config->cells[atom_i.cell_idx];

    double rho_val = 0.0; // initialize density for this atom
    double dF = 0.0;      // initialize gradient of embedding fn for this atom

    // Loop over pairs for this atom
    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairEAM &pair_ij = *(static_cast<PairEAM *>(pair_ij_ptr));  // tmp pair

      // Check that neighbor length lies in pair potential radius
      if (pair_ij.phi_knot != -1) {
        AtomEAM &atom_j = *(static_cast<AtomEAM *>(pair_ij.neigh));  // tmp atom

        // Compute phi(r_ij) and its gradient in one step
        double phigrad;
        double phival = 0.5 * phi_fns[pair_ij.phi_idx]->eval_comb(pair_ij.r, &phigrad);

        phigrad *= 0.5; // only half of the gradient/energy contributes to the force/energy since we are double counting

        cell.energy += phival;  // add in piece contributed by neighbor to energy

        Vect tmp_force = pair_ij.dist * phigrad;  // compute tmp force values
        atom_i.force += tmp_force;  // add in force on atom i from atom j
        atom_j.force -= tmp_force;  // subtract off force on atom j from atom i (Newton's law: action = -reaction)

        // Compute stress on cell
        tmp_force *= pair_ij.r;
        cell.stress -= pair_ij.dist & tmp_force;
      } // END IF STMNT: PAIR LIES INSIDE CUTOFF FOR PAIR POTENTIAL

      // Check that neighbor length lies in rho potential (density function) radius
      if (pair_ij.rho_knot != -1) {
        // Compute density and its gradient in one step
        rho_val += rho_fns[pair_ij.rho_idx]->eval_comb(pair_ij.r, &pair_ij.drho);
      } else {
        pair_ij.drho = 0.0;
      } // END IF STMNT: PAIR LIES INSIDE CUTOFF FOR RHO POTENTIAL
    } // END LOOP OVER PAIRS

    // Compute energy, gradient for embedding function F
    // Punish this potential for having rho lie outside of F
    if ( rho_val < F_fns[atom_i.F_idx]->get_min_rcut() ) {
      double rho_i = F_fns[atom_i.F_idx]->get_min_rcut();
      constraint_err[atom_i.cell_idx] += cell.weight * DUMMY_WEIGHT * 10. * (rho_i - rho_val) * (rho_i - rho_val);
      if (!embed_extrap_)
        rho_val = rho_i;  // set the density to the inner cutoff if we don't extrapolate embedding fn later
    } else if ( rho_val > F_fns[atom_i.F_idx]->get_max_rcut() ) {
      double rho_f = F_fns[atom_i.F_idx]->get_max_rcut();
      constraint_err[atom_i.cell_idx] += cell.weight * DUMMY_WEIGHT * 10. * (rho_val - rho_f) * (rho_val - rho_f);
      if (!embed_extrap_)
        rho_val = rho_f;  // set the density to the outer cutoff if we don't extrapolate embedding fn later
    }

    // Add energy contribution from embedding function and get gradient in one step
    cell.energy += F_fns[atom_i.F_idx]->eval_comb(rho_val, &dF);

    // Loop over pairs for this atom to compute EAM force
    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairEAM &pair_ij = *(static_cast<PairEAM *>(pair_ij_ptr));  // tmp pair
      AtomEAM &atom_j = *(static_cast<AtomEAM *>(pair_ij.neigh));  // tmp atom

      Vect tmp_force = pair_ij.dist * pair_ij.drho * dF;  // compute tmp force values
      atom_i.force += tmp_force;  // add in force on atom i from atom j
      atom_j.force -= tmp_force;  // subtract off force on atom j from atom i (Newton's law: action = -reaction)

      // Compute stress on cell
      tmp_force *= pair_ij.r;
      cell.stress -= pair_ij.dist & tmp_force;
    } // END 2nd LOOP OVER PAIRS
  } // END 1st LOOP OVER ATOMS

  accumulate_error(comm, error_vec, constraint_err);

  // Punishment for U'(n_mean) != 0
  for (int i=0; i<mmz->potlist->get_ntypes(); ++i) {
	  double rho_i = F_fns[i]->get_min_rcut();
	  double rho_f = F_fns[i]->get_max_rcut();
    double eam_error = DUMMY_WEIGHT * F_fns[i]->eval_grad(0.5 * (rho_i + rho_f));
    error_sum_ += eam_error * eam_error;
    if (error_vec && comm.is_root()) error_vec->push_back(eam_error);
  }

  ++ncalls_;  // keep track of the number of times this function is called

  return error_sum_;
}

/* ----------------------------------------------------------------------
   rescale potential
------------------------------------------------------------------------- */

int PotEAM::rescale(const Comm& comm, std::ostream *out, int flag)
{
  if (flag < 0) return 0;

  // Store embedding function
  std::vector<Basis *> F_fns;
  for (Basis*& fn : pot_fns_[2].fns)
    F_fns.push_back(fn);
  int F_size = F_fns.size();

  // Get atomic densities
  Vector<double> densities;
  compute_densities(comm, &densities);

  if (!comm.is_root()) return 0;

  // Store max and min densities for each embedding function
  std::vector<double> minrho(F_size,1e100), maxrho(F_size,-1e100);
  for (int i=0; i<mmz->config->total_natoms; ++i) {
    int typ = mmz->config->atoms[i]->typ;
    maxrho[typ] = std::max(densities[i], maxrho[typ]);
    minrho[typ] = std::min(densities[i], minrho[typ]);
  }

  // Loop through each alloy's embedding functions looking for global max/min potential idx
  double min_density = 1e100, max_density = -1e100;
  int min_F_idx = -1, max_F_idx = -1;
  for (int i=0; i<F_size; ++i) {
    if (maxrho[i]>max_density) {
      max_density = maxrho[i];
      max_F_idx = i;
    }
    if (minrho[i]<min_density) {
      min_density = minrho[i];
      min_F_idx = i;
    }
  }

  int sign = (max_density >= -min_density) ? 1 : -1;  // Determine the dominant side

  // Determine new left and right boundary, add 40 percent...
  // Loop through each column of F
  Vector<double> left(F_size), right(F_size);
  for (int i=0; i<F_size; ++i) {
    double rho_i = F_fns[i]->get_min_rcut();
    double rho_f = F_fns[i]->get_max_rcut();
    double length = rho_f - rho_i;
    left[i]  = minrho[i] - 0.01 * length; // left  boundary that just encompasses minrho by 1% of length
    right[i] = maxrho[i] + 0.01 * length; // right boundary that just encompasses maxrho by 1% of length

    if (flag || minrho[i] - rho_i < 0. ||       // is minrho/maxrho outside of range of F's axis
        minrho[i] - rho_i > 0.05 * length ||    // or is it too far inside the range of F's axis
        maxrho[i] - rho_f > 0. ||
        maxrho[i] - rho_f < -0.05 * length) flag = 1; // continue with scaling

    // Don't allow rescaling if bounds are too big
    if (right[i] - rho_f > 100.0 || left[i] - rho_i < -100.0) return 0;
  }

  double a = (sign==1) ? 1.0/right[max_F_idx] : 1.0/left[min_F_idx];  // scaling factor

  if (flag || std::abs(a)>1.05 || std::abs(a)<0.95) flag = 1; // scale if scaling factor increase/decrease pot by 5%

  if (!flag) return 0;

  // Rescale embedding fn with new x-axis boundaries
  for (int i=0; i<F_size; ++i) {
    F_fns[i]->readjust_x(left[i], right[i]);  // set new boundaries
    F_fns[i]->rescale_x(a); // rescale the x-coords by factor
  }

  // Scale density function to keep physics same
  std::vector<Basis *> rho_fns = pot_fns_[1].fns;
  for (Basis *rho_fn : rho_fns)
    *rho_fn *= a;

  // Output scaling factor to screen
  if (comm.is_root() && out)
    *out << "Embedding function scaling factor " << std::fixed << a << std::endl;

  // Reset 2nd derivs for all potentials
  for (PotFns& pot_fn : pot_fns_)
    pot_fn.refresh_basis();

  // Correct gauge: U'(n_mean)=0
  // phi_ij(r) -> phi_ij(r) + lambda_i*rho_j(r) + lambda_j*rho_i(r)
  // U_i(n_i) -> U_i(n_i) - lambda_i*n_i
  // So we need U_i'(n_i-mean) = lambda_i to satisfy first line

  // Define lambda for each atom type (or number of embedding fn's)
  std::vector<double> lambda(F_size);
  for (int i=0; i<F_size; ++i) {
    double rho_i = F_fns[i]->get_min_rcut();
    double rho_f = F_fns[i]->get_max_rcut();
    lambda[i] = F_fns[i]->eval_grad(0.5 * (rho_i + rho_f));
  }

  // phi_ij(r) -> phi_ij(r) + lambda_i*rho_j(r) + lambda_j*rho_i(r)
  std::vector<Basis *> phi_fns = pot_fns_[0].fns;
  int count = 0;
  for (int i=0; i<F_size; ++i) {
    for (int j=i; j<F_size; ++j) {
      Basis *tmp_basis;

      tmp_basis = rho_fns[i]->clone();
      *tmp_basis *= lambda[j];
      *phi_fns[count] += *tmp_basis;
      delete tmp_basis;

      tmp_basis = rho_fns[j]->clone();
      *tmp_basis *= lambda[i];
      *phi_fns[count] += *tmp_basis;
      delete tmp_basis;

      ++count;
    }
  }

  // U_i(n_i) -> U_i(n_i) - lambda_i*n_i
  for (int i=0; i<F_size; ++i)
    F_fns[i]->add_linear(-1.0 * lambda[i]);

  // Output each lambda to screen
  if (comm.is_root() && out)
    for (int i=0; i<F_size; ++i)
      *out << "lambda[" << i << "] = " << std::fixed << lambda[i] << std::endl;

  // Reset 2nd derivs for all potentials
  for (PotFns& pot_fn : pot_fns_)
    pot_fn.refresh_basis();

  return 1;
}

/* ----------------------------------------------------------------------
   compute atomic densities
------------------------------------------------------------------------- */

void PotEAM::compute_densities(const Comm& comm, Vector<double> *density_ptr)
{
  // Untrap procs if necessary
  if (is_trapped_) {
    int flag = 4;
    comm.bcast(&flag, 1, MPI_INT, comm.get_root());
  }

  // Initialize potential on all procs
  initialize_pot(comm);

  // Setup local variables for potential function
  std::vector<Basis *> rho_fns;
  for (Basis*& fn : pot_fns_[1].fns)
    rho_fns.push_back(fn);

  // Make list of densities for each atom
  int natoms = mmz->config->total_natoms;
  Vector<double> densities(natoms, 0.0);

  // Loop over all atoms in atomvec
  for (Atom*& atom_i_ptr : mmz->atomvec->atoms) {
    AtomEAM &atom_i = *(static_cast<AtomEAM *>(atom_i_ptr));  // tmp atom

    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairEAM &pair_ij = *(static_cast<PairEAM *>(pair_ij_ptr));  // tmp pair

      if (pair_ij.rho_knot != -1)  // pair distance is inside density function
        densities[atom_i.global_idx] += rho_fns[pair_ij.rho_idx]->eval(pair_ij.r);
    } // END LOOP OVER PAIRS
  } // END LOOP OVER ATOMS


  // Gather up densities from all procs
  std::vector<double> densities_final(natoms, 0.0);
  comm.reduce(&densities[0], &densities_final[0], natoms, MPI_DOUBLE, MPI_SUM, comm.get_root());

  if (comm.is_root()) density_ptr->swap(densities_final);
}

/* ----------------------------------------------------------------------
   write punishment/constraint errors (none for this pot type)
------------------------------------------------------------------------- */

void PotEAM::write_punish(const ErrorVec& e_vec, String filename, int num) const
{
  filename += ".punish";
  std::ostringstream oss;
  oss << filename;
  if (num >= 0)
  oss << "." << num;

  std::ofstream punish_fs (oss.str(), std::ios_base::out);
  punish_fs << "Limiting constraints (atomic density lies inside F boundaries)" << std::endl;
  punish_fs << std::setw(3) << std::right << "#"
            << std::setw(11) << std::right << "conf_w"
            << std::setw(15) << std::right << "punish error"
            << std::setw(15) << std::right << "punishment"
            << std::endl;

  // Error vector contains the following initial components
  // (1 energy, 6 stresses)*ncells = 7*ncells
  // (3 forces per atom)*natoms = 3*natoms
  int punish_idx = 7*mmz->config->ncells + 3*mmz->config->total_natoms;

  for (int c=0; c<mmz->config->ncells; ++c) {
    Cell &cell = *mmz->config->cells[c];
    punish_fs << std::setw(3) << std::right << c << std::fixed
              << std::setw(11) << std::right << std::setprecision(4) << cell.weight
              << std::setw(15) << std::right << std::setprecision(6) << e_vec[punish_idx]*e_vec[punish_idx]
              << std::setw(15) << std::right << std::setprecision(6) << e_vec[punish_idx]/cell.weight
              << std::endl;
    ++punish_idx;
  }

  punish_fs << std::endl << "Dummy constraints (F'(n_avg) = 0)" << std::endl;
  punish_fs << std::setw(7) << std::right << "element"
            << std::setw(15) << std::right << "punish error"
            << std::setw(15) << std::right << "F'(n_avg)"
            << std::endl;

  int ntypes = mmz->potlist->get_ntypes();
  int dummy_idx = punish_idx;
  for (int i=0; i<ntypes; ++i) {
    punish_fs << std::setw(7) << std::right << i << std::fixed
              << std::setw(15) << std::right << std::setprecision(6) << e_vec[dummy_idx]*e_vec[dummy_idx]
              << std::setw(15) << std::right << std::setprecision(6) << e_vec[dummy_idx]/DUMMY_WEIGHT
              << std::endl;
    ++dummy_idx;
  }

  punish_fs.close();

  return;
}

/* ----------------------------------------------------------------------
  write out extra data (densities for this pot)
------------------------------------------------------------------------- */

void PotEAM::write_extras(const Comm& comm, String filename, int num)
{
  write_densities(comm, filename, num);

  return;
}

/* ----------------------------------------------------------------------
  write out densities
------------------------------------------------------------------------- */

void PotEAM::write_densities(const Comm& comm, String filename, int num)
{
  // Get atomic densities
  std::vector<double> densities;
  compute_densities(comm, &densities);

  if (!comm.is_root()) return;

  filename += ".rho";
  std::ostringstream oss;
  oss << filename;
  if (num >= 0)
  oss << "." << num;

  std::ofstream rho_fs (oss.str(), std::ios_base::out);
  rho_fs << std::setw(6) << std::right << "conf:"
         << std::setw(5) << std::right << "atom"
         << std::setw(8) << std::right << "type"
         << std::setw(15) << std::right << "rho"
         << std::endl;

  int counter = -1;
  for (int c=0; c<mmz->config->ncells; ++c) {
    for (int i=0; i<mmz->config->cells[c]->natoms; ++i) {
      Atom &atom = *mmz->config->cells[c]->atoms[i];

      rho_fs << std::setw(5) << std::right << c << ":"
             << std::setw(5) << std::right << i
             << std::setw(8) << std::right << atom.typ
             << std::setw(15) << std::right << std::fixed << densities[++counter]
             << std::endl;
    }
  }

  rho_fs.close();

  return;
}
