#include "pot_pair_spline.h"
#include "atom_pair_spline.h"
#include "pair_pair_spline.h"
#include "../../config.h"
#include "../../cell.h"
#include "../../atom_vec.h"
#include "../../style_basis.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

PotPairSpline::PotPairSpline(class Meamzilla *mmz, int ntypes) : PotPair(mmz, ntypes)
{
  // ctor
}

/* ---------------------------------------------------------------------- */

PotPairSpline::~PotPairSpline()
{
  // dtor
}

/* ----------------------------------------------------------------------
   clone potential
------------------------------------------------------------------------- */

PotPairSpline* PotPairSpline::clone() const
{
  return new PotPairSpline(*this);
}

/* ----------------------------------------------------------------------
   compute error sum/vector
------------------------------------------------------------------------- */

double PotPairSpline::compute(const Comm& comm, ErrorVec *error_vec)
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
   compute error sum/vector
------------------------------------------------------------------------- */

template <class T>
double PotPairSpline::fast_compute(const Comm& comm, ErrorVec *error_vec)
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

  // 1st Loop over all atoms in atomvec
  for (Atom*& atom_i_ptr : mmz->atomvec->atoms) {
    // Make temporary atom and cell
    AtomPairSpline &atom_i = *(static_cast<AtomPairSpline*>(atom_i_ptr));
    Cell &cell = *mmz->config->cells[atom_i.cell_idx];

    // Loop over pairs for this atom
    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      // Make temporary pair and atom
      PairPairSpline &pair_ij = *(static_cast<PairPairSpline *>(pair_ij_ptr));
      AtomPairSpline &atom_j  = *(static_cast<AtomPairSpline *>(pair_ij.neigh));

      int self = (&atom_i == &atom_j) ? 1 : 0;

      // Check that neighbor length lies in pair potential radius
      if (pair_ij.phi_knot != -1) {
        // Compute phi(r_ij) and its gradient in one step
        double phigrad;
        double phival = phi_fns[pair_ij.phi_idx]->T::splint_comb(pair_ij.phi_knot, pair_ij.phi_shift, &phigrad);

        // Only half of the value and gradient contributes to the energy
        // and force if atom is self to avoid double counting
        if (self) {
          phival  *= 0.5;
          phigrad *= 0.5;
        }

        cell.energy += phival;  // add in piece contributed by neighbor to energy

        Vect tmp_force = pair_ij.dist * phigrad;  // compute tmp force values
        atom_i.force += tmp_force;  // add in force on atom i from atom j
        atom_j.force -= tmp_force;  // subtract off force on atom j from atom i (Newton's law: action = -reaction)

        // Compute stress on cell
        tmp_force *= pair_ij.r;
        cell.stress -= pair_ij.dist & tmp_force;
      } // END IF STMNT: PAIR LIES INSIDE CUTOFF FOR PAIR POTENTIAL
    } // END LOOP OVER PAIRS
  } // END 1st LOOP OVER ATOMS

  ErrorVec extra_error(0);
  accumulate_error(comm, error_vec, extra_error);

  ++ncalls_;  // keep track of the number of times this function is called

  return error_sum_;
}
