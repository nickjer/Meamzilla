#include "pot_jmeam2_spline.h"
#include "atom_jmeam2_spline.h"
#include "pair_jmeam2_spline.h"
#include "triplet_jmeam2_spline.h"
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

PotJMEAM2Spline::PotJMEAM2Spline(class Meamzilla *mmz, int ntypes) : PotJMEAM2(mmz, ntypes)
{
  // ctor
}

/* ---------------------------------------------------------------------- */

PotJMEAM2Spline::~PotJMEAM2Spline()
{
  // dtor
}

/* ----------------------------------------------------------------------
   clone potential
------------------------------------------------------------------------- */

PotJMEAM2Spline* PotJMEAM2Spline::clone() const
{
  return new PotJMEAM2Spline(*this);
}

/* ----------------------------------------------------------------------
   compute error sum/vector
------------------------------------------------------------------------- */

double PotJMEAM2Spline::compute(const Comm& comm, ErrorVec *error_vec)
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

void PotJMEAM2Spline::compute_densities(const Comm& comm, Vector<double> *density_ptr)
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
   fast compute error sum/vector
------------------------------------------------------------------------- */

template <class T>
double PotJMEAM2Spline::fast_compute(const Comm& comm, ErrorVec *error_vec)
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
  std::vector<T *> f_fns;
  for (Basis*& fn : pot_fns_[3].fns)
    f_fns.push_back(static_cast<T *>(fn));
  std::vector<T *> g_fns;
  for (Basis*& fn : pot_fns_[4].fns)
    g_fns.push_back(static_cast<T *>(fn));
  std::vector<T *> p_fns;
  for (Basis*& fn : pot_fns_[5].fns)
    p_fns.push_back(static_cast<T *>(fn));
  std::vector<T *> q_fns;
  for (Basis*& fn : pot_fns_[6].fns)
    q_fns.push_back(static_cast<T *>(fn));
  std::vector<T *> p2_fns;
  for (Basis*& fn : pot_fns_[7].fns)
    p2_fns.push_back(static_cast<T *>(fn));
  std::vector<T *> q2_fns;
  for (Basis*& fn : pot_fns_[8].fns)
    q2_fns.push_back(static_cast<T *>(fn));

  // Set up constraint error (error from density going out of bounds of embedding function)
  Vector<double> constraint_err(mmz->config->ncells,0.0);

  // Loop over all atoms in atomvec
  for (Atom*& atom_i_ptr : mmz->atomvec->atoms) {
    // Make temporary atom and cell
    AtomJMEAM2Spline &atom_i = *(static_cast<AtomJMEAM2Spline *>(atom_i_ptr));
    Cell &cell = *mmz->config->cells[atom_i.cell_idx];

    double rho_val = 0.0; // initialize density for this atom
    double dF = 0.0;      // initialize gradient of embedding fn for this atom

    // Loop over pairs for this atom
    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairJMEAM2Spline &pair_ij = *(static_cast<PairJMEAM2Spline *>(pair_ij_ptr));  // tmp pair

      // Check that neighbor length lies in pair potential radius
      if (pair_ij.phi_knot != -1) {
        AtomJMEAM2Spline &atom_j = *(static_cast<AtomJMEAM2Spline *>(pair_ij.neigh));  // tmp atom

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

      // Check that neighbor length lies in f- potential radius
      if (pair_ij.f_knot != -1) {
        pair_ij.f = f_fns[pair_ij.f_idx]->T::splint_comb(pair_ij.f_knot, pair_ij.f_shift, &pair_ij.df);
      } else {
        pair_ij.f = 0.0;
        pair_ij.df = 0.0;
      } // END IF STMNT: PAIR LIES INSIDE CUTOFF FOR f- POTENTIAL

      // Check that neighbor length lies in p- potential radius
      if (pair_ij.p_knot != -1) {
        pair_ij.p = p_fns[pair_ij.p_idx]->T::splint_comb(pair_ij.p_knot, pair_ij.p_shift, &pair_ij.dp);
      } else {
        pair_ij.p = 0.0;
        pair_ij.dp = 0.0;
      } // END IF STMNT: PAIR LIES INSIDE CUTOFF FOR p- POTENTIAL

      // Check that neighbor length lies in p- potential radius
      if (pair_ij.p2_knot != -1) {
        pair_ij.p2 = p2_fns[pair_ij.p2_idx]->T::splint_comb(pair_ij.p2_knot, pair_ij.p2_shift, &pair_ij.dp2);
      } else {
        pair_ij.p2 = 0.0;
        pair_ij.dp2 = 0.0;
      } // END IF STMNT: PAIR LIES INSIDE CUTOFF FOR p- POTENTIAL
    } // END LOOP OVER PAIRS

    // Loop over every angle formed by pairs called triplets
    for (Triplet*& triplet_ijk_ptr : atom_i.triplets) {
      TripletJMEAM2Spline &triplet_ijk = *(static_cast<TripletJMEAM2Spline *>(triplet_ijk_ptr));  // tmp triplet
      PairJMEAM2Spline &pair_ij = *(static_cast<PairJMEAM2Spline *>(triplet_ijk.pair_ij));  // tmp pairs
      PairJMEAM2Spline &pair_ik = *(static_cast<PairJMEAM2Spline *>(triplet_ijk.pair_ik));

      // The cos(theta) should always lie inside -1 ... 1
      // So store the g and g' without checking bounds
      triplet_ijk.g = g_fns[triplet_ijk.g_idx]->T::splint_comb(triplet_ijk.g_knot, triplet_ijk.g_shift, &triplet_ijk.dg);
      triplet_ijk.q = q_fns[triplet_ijk.q_idx]->T::splint_comb(triplet_ijk.q_knot, triplet_ijk.q_shift, &triplet_ijk.dq);
      triplet_ijk.q2 = q2_fns[triplet_ijk.q2_idx]->T::splint_comb(triplet_ijk.q2_knot, triplet_ijk.q2_shift, &triplet_ijk.dq2);

      // Sum up rho piece for atom i caused by j and k
      // f_ij * f_ik * g_ijk
      rho_val     += pair_ij.f * pair_ik.f * triplet_ijk.g;
      cell.energy += pair_ij.p * pair_ik.p * triplet_ijk.q;
      cell.energy += pair_ij.p2 * pair_ik.p2 * triplet_ijk.q2;
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
    cell.energy += F_fns[atom_i.F_idx]->T::splint_comb(rho_val, &dF);

    // Loop over pairs for this atom to compute EAM force
    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairJMEAM2Spline &pair_ij = *(static_cast<PairJMEAM2Spline *>(pair_ij_ptr));  // tmp pair
      AtomJMEAM2Spline &atom_j = *(static_cast<AtomJMEAM2Spline *>(pair_ij.neigh));  // tmp atom

      Vect tmp_force = pair_ij.dist * pair_ij.drho * dF;  // compute tmp force values
      atom_i.force += tmp_force;  // add in force on atom i from atom j
      atom_j.force -= tmp_force;  // subtract off force on atom j from atom i (Newton's law: action = -reaction)

      // Compute stress on cell
      tmp_force *= pair_ij.r;
      cell.stress -= pair_ij.dist & tmp_force;
    } // END 2nd LOOP OVER PAIRS

    // Loop over every angle formed by pairs called triplets
    for (Triplet*& triplet_ijk_ptr : atom_i.triplets) {
      TripletJMEAM2Spline &triplet_ijk = *(static_cast<TripletJMEAM2Spline *>(triplet_ijk_ptr));  // tmp triplet
      PairJMEAM2Spline &pair_ij = *(static_cast<PairJMEAM2Spline *>(triplet_ijk.pair_ij));  // tmp pairs
      PairJMEAM2Spline &pair_ik = *(static_cast<PairJMEAM2Spline *>(triplet_ijk.pair_ik));
      AtomJMEAM2Spline &atom_j = *(static_cast<AtomJMEAM2Spline *>(pair_ij.neigh));  // tmp atoms
      AtomJMEAM2Spline &atom_k = *(static_cast<AtomJMEAM2Spline *>(pair_ik.neigh));

      // Some tmp variables to clean up force fn below
      double dV3j = triplet_ijk.g  * pair_ij.df * pair_ik.f  * dF + triplet_ijk.q  * pair_ij.dp * pair_ik.p  + triplet_ijk.q2  * pair_ij.dp2 * pair_ik.p2;
      double dV3k = triplet_ijk.g  * pair_ij.f  * pair_ik.df * dF + triplet_ijk.q  * pair_ij.p  * pair_ik.dp + triplet_ijk.q2  * pair_ij.p2  * pair_ik.dp2;
      double V3   = triplet_ijk.dg * pair_ij.f  * pair_ik.f  * dF + triplet_ijk.dq * pair_ij.p  * pair_ik.p  + triplet_ijk.dq2 * pair_ij.p2  * pair_ik.p2;

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

  // Punishment for p-pot y-max magnitude not being 1.0
  double max_p_mag = std::abs(pot_fns_[5].get_max_y_mag());
  double p_pot_error = DUMMY_WEIGHT * 25. * (1.0 - max_p_mag) * (1.0 - max_p_mag);
  error_sum_ += p_pot_error * p_pot_error;
  if (error_vec && comm.is_root()) error_vec->push_back(p_pot_error);

  // Punishment for p2-pot y-max magnitude not being 1.0
  double max_p2_mag = std::abs(pot_fns_[7].get_max_y_mag());
  double p2_pot_error = DUMMY_WEIGHT * 25. * (1.0 - max_p2_mag) * (1.0 - max_p2_mag);
  error_sum_ += p2_pot_error * p2_pot_error;
  if (error_vec && comm.is_root()) error_vec->push_back(p2_pot_error);

  ++ncalls_;  // keep track of the number of times this function is called

  return error_sum_;
}

/* ----------------------------------------------------------------------
   compute atomic densities
------------------------------------------------------------------------- */

template <class T>
void PotJMEAM2Spline::fast_compute_densities(const Comm& comm, Vector<double> *density_ptr)
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
  std::vector<T *> f_fns;
  for (Basis*& fn : pot_fns_[3].fns)
    f_fns.push_back(static_cast<T *>(fn));
  std::vector<T *> g_fns;
  for (Basis*& fn : pot_fns_[4].fns)
    g_fns.push_back(static_cast<T *>(fn));

  // Make list of densities for each atom
  int natoms = mmz->config->total_natoms;
  Vector<double> densities(natoms, 0.0);

  // Loop over all atoms in atomvec
  for (Atom*& atom_i_ptr : mmz->atomvec->atoms) {
    AtomJMEAM2Spline &atom_i = *(static_cast<AtomJMEAM2Spline *>(atom_i_ptr));  // tmp atom

    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairJMEAM2Spline &pair_ij = *(static_cast<PairJMEAM2Spline *>(pair_ij_ptr));  // tmp pair

      if (pair_ij.rho_knot != -1)  // pair distance is inside density function
        densities[atom_i.global_idx] += rho_fns[pair_ij.rho_idx]->T::splint(pair_ij.rho_knot, pair_ij.rho_shift);

      if (pair_ij.f_knot != -1)  // radial distance inside f-potential
        pair_ij.f = f_fns[pair_ij.f_idx]->T::splint(pair_ij.f_knot, pair_ij.f_shift);
      else
        pair_ij.f = 0.0;
    } // END LOOP OVER PAIRS

    for (Triplet*& triplet_ijk_ptr : atom_i.triplets) {
      TripletJMEAM2Spline &triplet_ijk = *(static_cast<TripletJMEAM2Spline *>(triplet_ijk_ptr));  // tmp triplet
      PairJMEAM2Spline &pair_ij = *(static_cast<PairJMEAM2Spline *>(triplet_ijk.pair_ij));  // tmp pairs
      PairJMEAM2Spline &pair_ik = *(static_cast<PairJMEAM2Spline *>(triplet_ijk.pair_ik));

      double g_val = g_fns[triplet_ijk.g_idx]->T::splint(triplet_ijk.g_knot, triplet_ijk.g_shift);

      densities[atom_i.global_idx] += pair_ij.f * pair_ik.f * g_val;
    } // END LOOP OVER TRIPLETS
  } // END LOOP OVER ATOMS


  // Gather up densities from all procs
  std::vector<double> densities_final(natoms, 0.0);
  comm.reduce(&densities[0], &densities_final[0], natoms, MPI_DOUBLE, MPI_SUM, comm.get_root());

  if (comm.is_root()) density_ptr->swap(densities_final);
}

