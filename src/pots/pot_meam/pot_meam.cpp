#include "pot_meam.h"
#include "atom_meam.h"
#include "pair_meam.h"
#include "triplet_meam.h"
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

PotMEAM::PotMEAM(class Meamzilla *mmz, int ntypes) : PotEAM(mmz, ntypes)
{
  // Double count pairs
  dbl_cnt_pair_ = 1;

  // Resize list of potentials
  npot_fns_ = 5;
  pot_fns_.resize(npot_fns_);

  // f- potential has n*(n+1)/2 potentials and is a radial potential
  // For types A and B, you have f_AA, f_AB, f_BB
  pot_fns_[3].setup_pot_fns(ntypes, PotFns::AlloyType::Pair_ij, PotFns::FnType::Radial);

  // g- potential has n potentials (depending on atom i, the original atom) and is NOT a radial potential
  // For types A and B, you have g_A, g_B
  pot_fns_[4].setup_pot_fns(ntypes, PotFns::AlloyType::Atom_i, PotFns::FnType::Undefined);
}

/* ---------------------------------------------------------------------- */

PotMEAM::~PotMEAM()
{
  // dtor
}

/* ----------------------------------------------------------------------
   clone potential
------------------------------------------------------------------------- */

PotMEAM* PotMEAM::clone() const
{
  return new PotMEAM(*this);
}

/* ----------------------------------------------------------------------
   check that a potential is similar enough to this
   potential for pairs/triplets to have identical info
------------------------------------------------------------------------- */

void PotMEAM::check_if_similar(Potential *pot_ptr) const
{
  PotMEAM& pot = *static_cast<PotMEAM *>(pot_ptr);

  // Compare phi
  pot_fns_[0].check_if_similar(pot.pot_fns_[0]);  // error will be thrown if not similar

  // Compare rho
  pot_fns_[1].check_if_similar(pot.pot_fns_[1]);  // error will be thrown if not similar

  // Do not compare embedding fn

  // Compare f
  pot_fns_[3].check_if_similar(pot.pot_fns_[3]);  // error will be thrown if not similar

  // Compare g
  pot_fns_[4].check_if_similar(pot.pot_fns_[4]);  // error will be thrown if not similar

  return;
}

/* ----------------------------------------------------------------------
   compute error sum/vector
------------------------------------------------------------------------- */

double PotMEAM::compute(const Comm& comm, ErrorVec *error_vec)
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
  std::vector<Basis *> f_fns;
  for (Basis*& fn : pot_fns_[3].fns)
    f_fns.push_back(fn);
  std::vector<Basis *> g_fns;
  for (Basis*& fn : pot_fns_[4].fns)
    g_fns.push_back(fn);

  // Set up constraint error (error from density going out of bounds of embedding function)
  Vector<double> constraint_err(mmz->config->ncells,0.0);

  // Loop over all atoms in atomvec
  for (Atom*& atom_i_ptr : mmz->atomvec->atoms) {
    // Make temporary atom and cell
    AtomMEAM &atom_i = *(static_cast<AtomMEAM *>(atom_i_ptr));
    Cell &cell = *mmz->config->cells[atom_i.cell_idx];

    double rho_val = 0.0; // initialize density for this atom
    double dF = 0.0;      // initialize gradient of embedding fn for this atom

    // Loop over pairs for this atom
    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairMEAM &pair_ij = *(static_cast<PairMEAM *>(pair_ij_ptr));  // tmp pair

      // Check that neighbor length lies in pair potential radius
      if (pair_ij.phi_knot != -1) {
        AtomMEAM &atom_j = *(static_cast<AtomMEAM *>(pair_ij.neigh));  // tmp atom

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

      // Check that neighbor length lies in f- potential radius
      if (pair_ij.f_knot != -1) {
        pair_ij.f = f_fns[pair_ij.f_idx]->eval_comb(pair_ij.r, &pair_ij.df);
      } else {
        pair_ij.f = 0.0;
        pair_ij.df = 0.0;
      } // END IF STMNT: PAIR LIES INSIDE CUTOFF FOR f- POTENTIAL
    } // END LOOP OVER PAIRS

    // Loop over every angle formed by pairs called triplets
    for (Triplet*& triplet_ijk_ptr : atom_i.triplets) {
      TripletMEAM &triplet_ijk = *(static_cast<TripletMEAM *>(triplet_ijk_ptr));  // tmp triplet
      PairMEAM &pair_ij = *(static_cast<PairMEAM *>(triplet_ijk.pair_ij));  // tmp pairs
      PairMEAM &pair_ik = *(static_cast<PairMEAM *>(triplet_ijk.pair_ik));

      // The cos(theta) should always lie inside -1 ... 1
      // So store the g and g' without checking bounds
      triplet_ijk.g = g_fns[triplet_ijk.g_idx]->eval_comb(triplet_ijk.cos, &triplet_ijk.dg);

      // Sum up rho piece for atom i caused by j and k
      // f_ij * f_ik * g_ijk
      rho_val += pair_ij.f * pair_ik.f * triplet_ijk.g;
    } // END LOOP OVER TRIPLETS

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
      PairMEAM &pair_ij = *(static_cast<PairMEAM *>(pair_ij_ptr));  // tmp pair
      AtomMEAM &atom_j = *(static_cast<AtomMEAM *>(pair_ij.neigh));  // tmp atom

      Vect tmp_force = pair_ij.dist * pair_ij.drho * dF;  // compute tmp force values
      atom_i.force += tmp_force;  // add in force on atom i from atom j
      atom_j.force -= tmp_force;  // subtract off force on atom j from atom i (Newton's law: action = -reaction)

      // Compute stress on cell
      tmp_force *= pair_ij.r;
      cell.stress -= pair_ij.dist & tmp_force;
    } // END 2nd LOOP OVER PAIRS

    // Loop over every angle formed by pairs called triplets
    for (Triplet*& triplet_ijk_ptr : atom_i.triplets) {
      TripletMEAM &triplet_ijk = *(static_cast<TripletMEAM *>(triplet_ijk_ptr));  // tmp triplet
      PairMEAM &pair_ij = *(static_cast<PairMEAM *>(triplet_ijk.pair_ij));  // tmp pairs
      PairMEAM &pair_ik = *(static_cast<PairMEAM *>(triplet_ijk.pair_ik));
      AtomMEAM &atom_j = *(static_cast<AtomMEAM *>(pair_ij.neigh));  // tmp atoms
      AtomMEAM &atom_k = *(static_cast<AtomMEAM *>(pair_ik.neigh));

      // Some tmp variables to clean up force fn below
      double dV3j = triplet_ijk.g  * pair_ij.df * pair_ik.f  * dF;
      double dV3k = triplet_ijk.g  * pair_ij.f  * pair_ik.df * dF;
      double V3   = triplet_ijk.dg * pair_ij.f  * pair_ik.f  * dF;

      double vlj  = V3 * pair_ij.invr;
      double vlk  = V3 * pair_ik.invr;
      double vv3j = dV3j - vlj * triplet_ijk.cos;
      double vv3k = dV3k - vlk * triplet_ijk.cos;

      Vect dfj = pair_ij.dist * vv3j + pair_ik.dist * vlj;
      Vect dfk = pair_ik.dist * vv3k + pair_ij.dist * vlk;

      atom_i.force += dfj + dfk;  // force on atom i from j and k
      atom_j.force -= dfj;  // reaction force on atom j from i and k
      atom_k.force -= dfk;  // reaction force on atom k from i and j

      // Compute stress on cell
      dfj *= pair_ij.r;
      dfk *= pair_ik.r;
      cell.stress -= pair_ij.dist & dfj;
      cell.stress -= pair_ik.dist & dfk;
    } // END LOOP OVER TRIPLETS
  } // END 1st LOOP OVER ATOMS

  accumulate_error(comm, error_vec, constraint_err);

  // Punishment for f-pot y-max magnitude not being 1.0
  double max_f_mag = std::abs(pot_fns_[3].get_max_y_mag());
  double f_pot_error = DUMMY_WEIGHT * 25. * (1.0 - max_f_mag) * (1.0 - max_f_mag);
  error_sum_ += f_pot_error * f_pot_error;
  if (error_vec && comm.is_root()) error_vec->push_back(f_pot_error);

  ++ncalls_;  // keep track of the number of times this function is called

  return error_sum_;
}

/* ----------------------------------------------------------------------
   rescale potential
------------------------------------------------------------------------- */

int PotMEAM::rescale(const Comm& comm, std::ostream *out, int flag)
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

  // Rescale MEAM potentials: f*f*g
  int flag_meam = rescale_3body(comm, out, flag);
  if (flag_meam && !flag) return 1;

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

  // In MEAM you have a*f*f*g, where 'a' scale factor is
  // carried over from scaling the total density
  // We multiply this 'a' to the g potential
  std::vector<Basis *> g_fns = pot_fns_[4].fns;
  for (Basis *g_fn : g_fns)
    *g_fn *= a;

  // Output scaling factor to screen
  if (comm.is_root() && out)
    *out << "Embedding function scaling factor " << std::fixed << a << std::endl;

  return 1;
}

/* ----------------------------------------------------------------------
   compute atomic densities
------------------------------------------------------------------------- */

void PotMEAM::compute_densities(const Comm& comm, Vector<double> *density_ptr)
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
  std::vector<Basis *> f_fns;
  for (Basis*& fn : pot_fns_[3].fns)
    f_fns.push_back(fn);
  std::vector<Basis *> g_fns;
  for (Basis*& fn : pot_fns_[4].fns)
    g_fns.push_back(fn);

  // Make list of densities for each atom
  int natoms = mmz->config->total_natoms;
  Vector<double> densities(natoms, 0.0);

  // Loop over all atoms in atomvec
  for (Atom*& atom_i_ptr : mmz->atomvec->atoms) {
    AtomMEAM &atom_i = *(static_cast<AtomMEAM *>(atom_i_ptr));  // tmp atom

    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairMEAM &pair_ij = *(static_cast<PairMEAM *>(pair_ij_ptr));  // tmp pair

      if (pair_ij.rho_knot != -1)  // pair distance is inside density function
        densities[atom_i.global_idx] += rho_fns[pair_ij.rho_idx]->eval(pair_ij.r);

      if (pair_ij.f_knot != -1)  // radial distance inside f-potential
        pair_ij.f = f_fns[pair_ij.f_idx]->eval(pair_ij.r);
      else
        pair_ij.f = 0.0;
    } // END LOOP OVER PAIRS

    for (Triplet*& triplet_ijk_ptr : atom_i.triplets) {
      TripletMEAM &triplet_ijk = *(static_cast<TripletMEAM *>(triplet_ijk_ptr));  // tmp triplet
      PairMEAM &pair_ij = *(static_cast<PairMEAM *>(triplet_ijk.pair_ij));  // tmp pairs
      PairMEAM &pair_ik = *(static_cast<PairMEAM *>(triplet_ijk.pair_ik));

      double g_val = g_fns[triplet_ijk.g_idx]->eval(triplet_ijk.cos);

      densities[atom_i.global_idx] += pair_ij.f * pair_ik.f * g_val;
    } // END LOOP OVER TRIPLETS
  } // END LOOP OVER ATOMS


  // Gather up densities from all procs
  std::vector<double> densities_final(natoms, 0.0);
  comm.reduce(&densities[0], &densities_final[0], natoms, MPI_DOUBLE, MPI_SUM, comm.get_root());

  if (comm.is_root()) density_ptr->swap(densities_final);
}

/* ----------------------------------------------------------------------
   write punishment/constraint errors (none for this pot type)
------------------------------------------------------------------------- */

void PotMEAM::write_punish(const ErrorVec& e_vec, String filename, int num) const
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

  punish_fs << std::endl << "Dummy constraints (|f_max(r)| = 1)" << std::endl;
  punish_fs << std::setw(15) << std::right << "punish error"
            << std::setw(15) << std::right << "|f_max(r)|"
            << std::endl;

  double max_f_mag = std::abs(pot_fns_[3].get_max_y_mag());

  punish_fs << std::fixed
            << std::setw(15) << std::right << std::setprecision(6) << e_vec[punish_idx]*e_vec[punish_idx]
            << std::setw(15) << std::right << std::setprecision(6) << std::abs(max_f_mag)
            << std::endl;
  ++punish_idx;

  punish_fs.close();

  return;
}

/* ----------------------------------------------------------------------
  rescale MEAM part of the potential: f*f*g
------------------------------------------------------------------------- */

int PotMEAM::rescale_3body(const Comm& comm, std::ostream *out, int flag)
{
  if (!flag) return 0; // Don't run rescaling if flag is 0

  // Get scaling factor
  double max_f_mag = pot_fns_[3].get_max_y_mag();
  double b = 1.0/max_f_mag;

  // Scale beta first: f*f*g = (b*f)(b*f)(g/b^2)
  // Scale f-pot first
  std::vector<Basis *> f_fns = pot_fns_[3].fns;
  for (Basis *f_fn : f_fns)
    *f_fn *= b;

  // We scale g by 1/b^2
  std::vector<Basis *> g_fns = pot_fns_[4].fns;
  for (Basis *g_fn : g_fns)
    *g_fn /= b*b;

  // Output scaling factor to screen
  if (comm.is_root() && out)
    *out << "MEAM potential scaling factor " << std::fixed << b << std::endl;

  return 1;
}
