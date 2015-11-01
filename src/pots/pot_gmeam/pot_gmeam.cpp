/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "pot_gmeam.h"
#include "../pot_meam/atom_meam.h"
#include "../pot_meam/pair_meam.h"
#include "../pot_meam/triplet_meam.h"
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

PotGMEAM::PotGMEAM(class Meamzilla *mmz, int ntypes) : PotMEAM(mmz, ntypes)
{
  // g- potential has n*n*(n+1)/2 potentials (depending on atom_i and the pair _j & _k) and is NOT a radial potential
  // For types A and B, you have g_AAA, g_AAB, g_ABB, g_BAA, g_BAB, g_BBB
  pot_fns_[4].setup_pot_fns(ntypes, PotFns::AlloyType::Triplet_i_jk, PotFns::FnType::Undefined);
}

/* ---------------------------------------------------------------------- */

PotGMEAM::~PotGMEAM()
{
  // dtor
}

/* ----------------------------------------------------------------------
   clone potential
------------------------------------------------------------------------- */

PotGMEAM* PotGMEAM::clone() const
{
  return new PotGMEAM(*this);
}

/* ----------------------------------------------------------------------
   compute error sum/vector
------------------------------------------------------------------------- */

double PotGMEAM::compute(const Comm& comm, ErrorVec *error_vec)
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
      rho_val     += pair_ij.f * pair_ik.f * triplet_ijk.g;
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
  for (Basis*& fn : pot_fns_[3].fns) {
    double max_f_mag = std::abs(fn->get_max_y_mag());
    double f_pot_error = DUMMY_WEIGHT * 25. * (1.0 - max_f_mag) * (1.0 - max_f_mag);
    error_sum_ += f_pot_error * f_pot_error;
    if (error_vec && comm.is_root()) error_vec->push_back(f_pot_error);
  }

  ++ncalls_;  // keep track of the number of times this function is called

  return error_sum_;
}

/* ----------------------------------------------------------------------
   write punishment/constraint errors (none for this pot type)
------------------------------------------------------------------------- */

void PotGMEAM::write_punish(const ErrorVec& e_vec, String filename, int num) const
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

  for (Basis* fn : pot_fns_[3].fns) {
    double max_f_mag = std::abs(fn->get_max_y_mag());

    punish_fs << std::fixed
              << std::setw(15) << std::right << std::setprecision(6) << e_vec[punish_idx]*e_vec[punish_idx]
              << std::setw(15) << std::right << std::setprecision(6) << std::abs(max_f_mag)
              << std::endl;
    ++punish_idx;
  }

  punish_fs.close();

  return;
}

/* ----------------------------------------------------------------------
  rescale MEAM part of the potential: f*f*g
------------------------------------------------------------------------- */

int PotGMEAM::rescale_3body(const Comm& comm, std::ostream *out, int flag)
{
  if (!flag) return 0; // Don't run rescaling if flag is 0

  int ntypes = mmz->potlist->get_ntypes();

  int f_idx = 0;
  for (Basis*& f_fn : pot_fns_[3].fns) {
    double max_f_mag = f_fn->get_max_y_mag();
    double b = 1.0/max_f_mag;

    // Scale f-pot
    *f_fn *= b;

    // Scale g-pot
    std::vector<Basis *> g_fns = pot_fns_[4].fns;
    for (int i=0; i<ntypes; ++i) {
      for (int j=0; j<ntypes; ++j) {
        for (int k=j; k<ntypes; ++k) {
          int ij_idx = pot_fns_[3].get_2body_alloy_idx(i, j);
          int ik_idx = pot_fns_[3].get_2body_alloy_idx(i, k);
          int ijk_idx = pot_fns_[4].get_3body_alloy_idx(i, j, k);

          if (f_idx == ij_idx)
            *g_fns[ijk_idx] /= b;
          if (f_idx == ik_idx)
            *g_fns[ijk_idx] /= b;
        }
      }
    }

    // Output scaling factor to screen
    if (comm.is_root() && out)
      *out << "MEAM potential scaling factor (b_" << f_idx << ") " << std::fixed << b << std::endl;

    ++f_idx;
  }

  return 1;
}
