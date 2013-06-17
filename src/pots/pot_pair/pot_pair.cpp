#include "pot_pair.h"
#include "atom_pair.h"
#include "pair_pair.h"
#include "../../config.h"
#include "../../cell.h"
#include "../../atom_vec.h"
#include "../../error.h"
#include "../../style_basis.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

PotPair::PotPair(class Meamzilla *mmz, int ntypes) : Potential(mmz, ntypes)
{
  // Double count pairs
  dbl_cnt_pair_ = 0;

  // Resize list of potentials
  npot_fns_ = 1;
  pot_fns_.resize(npot_fns_);

  // Pair potential has n*(n+1)/2 potentials and is a radial potential
  // For types A and B, you have phi_AA, phi_AB, phi_BB
  pot_fns_[0].setup_pot_fns(ntypes, PotFns::AlloyType::Pair_ij, PotFns::FnType::Radial);
}

/* ---------------------------------------------------------------------- */

PotPair::~PotPair()
{
  // dtor
}

/* ----------------------------------------------------------------------
   clone potential
------------------------------------------------------------------------- */

PotPair* PotPair::clone() const
{
  return new PotPair(*this);
}

/* ----------------------------------------------------------------------
   check that a potential is similar enough to this
   potential for pairs/triplets to have identical info
------------------------------------------------------------------------- */

void PotPair::check_if_similar(Potential *pot_ptr) const
{
  PotPair& pot = *static_cast<PotPair *>(pot_ptr);

  // Compare phi
  pot_fns_[0].check_if_similar(pot.pot_fns_[0]);  // error will be thrown if not similar

  return;
}

/* ----------------------------------------------------------------------
   compute error sum/vector
------------------------------------------------------------------------- */

double PotPair::compute(const Comm& comm, ErrorVec *error_vec)
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

  // 1st Loop over all atoms in atomvec
  for (Atom*& atom_i_ptr : mmz->atomvec->atoms) {
    // Make temporary atom and cell
    AtomPair &atom_i = *(static_cast<AtomPair*>(atom_i_ptr));
    Cell &cell = *mmz->config->cells[atom_i.cell_idx];

    // Loop over pairs for this atom
    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      // Make temporary pair and atom
      PairPair &pair_ij = *(static_cast<PairPair *>(pair_ij_ptr));
      AtomPair &atom_j  = *(static_cast<AtomPair *>(pair_ij.neigh));

      int self = (&atom_i == &atom_j) ? 1 : 0;

      // Check that neighbor length lies in pair potential radius
      if (pair_ij.phi_knot != -1) {
        // Compute phi(r_ij) and its gradient in one step
        double phigrad;
        double phival = phi_fns[pair_ij.phi_idx]->eval_comb(pair_ij.r, &phigrad);

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
