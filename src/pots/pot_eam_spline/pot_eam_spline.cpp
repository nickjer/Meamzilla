#include "pot_eam_spline.h"
#include "atom_eam_spline.h"
#include "pair_eam_spline.h"
#include "../../input.h"
#include "../../config.h"
#include "../../cell.h"
#include "../../atom_vec.h"
#include "../../error.h"
#include "../../pot_list.h"
#include "../../style_basis.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

PotEAMSpline::PotEAMSpline(class Meamzilla *mmz, int ntypes) : PotEAM(mmz, ntypes)
{
  // ctor
}

/* ---------------------------------------------------------------------- */

PotEAMSpline::~PotEAMSpline()
{
  // dtor
}

/* ----------------------------------------------------------------------
   compute error sum/vector
------------------------------------------------------------------------- */

double PotEAMSpline::compute(const Comm& comm, ErrorVec *error_vec)
{
  String common_basis_type = pot_fns_[0].get_common_basis_type();

  for (PotFns& pot_fn : pot_fns_)
    if (common_basis_type != pot_fn.get_common_basis_type())
      common_basis_type = "undefined";

  if (common_basis_type == "undefined") fast_compute<Spline>(comm, error_vec);  // default
#define BASIS_CLASS
#define BasisStyle(key,Class) \
  else if (common_basis_type == #key) fast_compute<Class>(comm, error_vec); // specific basis
#include "../../style_basis.h"
#undef BasisStyle
#undef BASIS_CLASS

  return error_sum_;
}

/* ----------------------------------------------------------------------
   compute atomic densities
------------------------------------------------------------------------- */

void PotEAMSpline::compute_densities(const Comm& comm, Vector<double> *density_ptr)
{
  String common_basis_type = pot_fns_[0].get_common_basis_type();

  for (PotFns& pot_fn : pot_fns_)
    if (common_basis_type != pot_fn.get_common_basis_type())
      common_basis_type = "undefined";

  if (common_basis_type == "undefined") fast_compute_densities<Spline>(comm, density_ptr);  // default
#define BASIS_CLASS
#define BasisStyle(key,Class) \
  else if (common_basis_type == #key) fast_compute_densities<Class>(comm, density_ptr); // specific basis
#include "../../style_basis.h"
#undef BasisStyle
#undef BASIS_CLASS

  return;
}

/* ----------------------------------------------------------------------
   clone potential
------------------------------------------------------------------------- */

PotEAMSpline* PotEAMSpline::clone() const
{
  return new PotEAMSpline(*this);
}

/* ----------------------------------------------------------------------
   compute error sum/vector
------------------------------------------------------------------------- */

template <class T>
double PotEAMSpline::fast_compute(const Comm& comm, ErrorVec *error_vec)
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
  std::vector<T *> phi_fns;
  for (Basis*& fn : pot_fns_[0].fns)
    phi_fns.push_back(static_cast<T *>(fn));
  std::vector<T *> rho_fns;
  for (Basis*& fn : pot_fns_[1].fns)
    rho_fns.push_back(static_cast<T *>(fn));
  std::vector<T *> F_fns;
  for (Basis*& fn : pot_fns_[2].fns)
    F_fns.push_back(static_cast<T *>(fn));

  // Set up constraint error (error from density going out of bounds of embedding function)
  Vector<double> constraint_err(mmz->config->ncells,0.0);

  // Loop over all atoms in atomvec
  for (Atom*& atom_i_ptr : mmz->atomvec->atoms) {
    // Make temporary atom and cell
    AtomEAMSpline &atom_i = *(static_cast<AtomEAMSpline *>(atom_i_ptr));
    Cell &cell = *mmz->config->cells[atom_i.cell_idx];

    double rho_val = 0.0; // initialize density for this atom
    double dF = 0.0;      // initialize gradient of embedding fn for this atom

    // Loop over pairs for this atom
    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairEAMSpline &pair_ij = *(static_cast<PairEAMSpline *>(pair_ij_ptr));  // tmp pair

      // Check that neighbor length lies in pair potential radius
      if (pair_ij.phi_knot != -1) {
        AtomEAMSpline &atom_j = *(static_cast<AtomEAMSpline *>(pair_ij.neigh));  // tmp atom

        // Compute phi(r_ij) and its gradient in one step
        double phigrad;
        double phival = 0.5 * phi_fns[pair_ij.phi_idx]->T::splint_comb(pair_ij.phi_knot, pair_ij.phi_shift, &phigrad);

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
        rho_val += rho_fns[pair_ij.rho_idx]->T::splint_comb(pair_ij.rho_knot, pair_ij.rho_shift, &pair_ij.drho);
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
    cell.energy += F_fns[atom_i.F_idx]->T::splint_comb(rho_val, &dF);

    // Loop over pairs for this atom to compute EAM force
    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairEAMSpline &pair_ij = *(static_cast<PairEAMSpline *>(pair_ij_ptr));  // tmp pair
      AtomEAMSpline &atom_j = *(static_cast<AtomEAMSpline *>(pair_ij.neigh));  // tmp atom

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
    double eam_error = DUMMY_WEIGHT * F_fns[i]->T::splint_grad(0.5 * (rho_i + rho_f));
    error_sum_ += eam_error * eam_error;
    if (error_vec && comm.is_root()) error_vec->push_back(eam_error);
  }

  ++ncalls_;  // keep track of the number of times this function is called

  return error_sum_;
}

/* ----------------------------------------------------------------------
   compute atomic densities
------------------------------------------------------------------------- */

template <class T>
void PotEAMSpline::fast_compute_densities(const Comm& comm, Vector<double> *density_ptr)
{
  // Untrap procs if necessary
  if (is_trapped_) {
    int flag = 4;
    comm.bcast(&flag, 1, MPI_INT, comm.get_root());
  }

  // Initialize potential on all procs
  initialize_pot(comm);

  // Setup local variables for potential function
  std::vector<T *> rho_fns;
  for (Basis*& fn : pot_fns_[1].fns)
    rho_fns.push_back(static_cast<T *>(fn));

  // Make list of densities for each atom
  int natoms = mmz->config->total_natoms;
  Vector<double> densities(natoms, 0.0);

  // Loop over all atoms in atomvec
  for (Atom*& atom_i_ptr : mmz->atomvec->atoms) {
    AtomEAMSpline &atom_i = *(static_cast<AtomEAMSpline *>(atom_i_ptr));  // tmp atom

    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairEAMSpline &pair_ij = *(static_cast<PairEAMSpline *>(pair_ij_ptr));  // tmp pair

      if (pair_ij.rho_knot != -1)  // pair distance is inside density function
        densities[atom_i.global_idx] += rho_fns[pair_ij.rho_idx]->T::splint(pair_ij.rho_knot, pair_ij.rho_shift);
    } // END LOOP OVER PAIRS
  } // END LOOP OVER ATOMS


  // Gather up densities from all procs
  std::vector<double> densities_final(natoms, 0.0);
  comm.reduce(&densities[0], &densities_final[0], natoms, MPI_DOUBLE, MPI_SUM, comm.get_root());

  if (comm.is_root()) density_ptr->swap(densities_final);
}
