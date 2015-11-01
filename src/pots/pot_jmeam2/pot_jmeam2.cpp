/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "pot_jmeam2.h"
#include "atom_jmeam2.h"
#include "pair_jmeam2.h"
#include "triplet_jmeam2.h"
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

PotJMEAM2::PotJMEAM2(class Meamzilla *mmz, int ntypes) : PotJMEAM(mmz, ntypes)
{
  // Double count pairs
  dbl_cnt_pair_ = 1;

  // Resize list of potentials
  npot_fns_ = 9;
  pot_fns_.resize(npot_fns_);

  // p- potential has n*(n+1)/2 potentials and is a radial potential
  // For types A and B, you have p_AA, p_AB, p_BB
  pot_fns_[7].setup_pot_fns(ntypes, PotFns::AlloyType::Pair_ij, PotFns::FnType::Radial);

  // q- potential has n potentials (depending on atom i, the original atom) and is NOT a radial potential
  // For types A and B, you have q_A, q_B
  pot_fns_[8].setup_pot_fns(ntypes, PotFns::AlloyType::Triplet_i, PotFns::FnType::Undefined);
}

/* ---------------------------------------------------------------------- */

PotJMEAM2::~PotJMEAM2()
{
  // dtor
}

/* ----------------------------------------------------------------------
   clone potential
------------------------------------------------------------------------- */

PotJMEAM2* PotJMEAM2::clone() const
{
  return new PotJMEAM2(*this);
}

/* ----------------------------------------------------------------------
   check that a potential is similar enough to this
   potential for pairs/triplets to have identical info
------------------------------------------------------------------------- */

void PotJMEAM2::check_if_similar(Potential *pot_ptr) const
{
  PotJMEAM2& pot = *static_cast<PotJMEAM2 *>(pot_ptr);

  // Compare phi
  pot_fns_[0].check_if_similar(pot.pot_fns_[0]);  // error will be thrown if not similar

  // Compare rho
  pot_fns_[1].check_if_similar(pot.pot_fns_[1]);  // error will be thrown if not similar

  // Do not compare embedding fn

  // Compare f
  pot_fns_[3].check_if_similar(pot.pot_fns_[3]);  // error will be thrown if not similar

  // Compare g
  pot_fns_[4].check_if_similar(pot.pot_fns_[4]);  // error will be thrown if not similar

  // Compare p
  pot_fns_[5].check_if_similar(pot.pot_fns_[5]);  // error will be thrown if not similar

  // Compare q
  pot_fns_[6].check_if_similar(pot.pot_fns_[6]);  // error will be thrown if not similar

  // Compare p2
  pot_fns_[7].check_if_similar(pot.pot_fns_[7]);  // error will be thrown if not similar

  // Compare q2
  pot_fns_[8].check_if_similar(pot.pot_fns_[8]);  // error will be thrown if not similar

  return;
}

/* ----------------------------------------------------------------------
   compute error sum/vector
------------------------------------------------------------------------- */

double PotJMEAM2::compute(const Comm& comm, ErrorVec *error_vec)
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
  std::vector<Basis *> p_fns;
  for (Basis*& fn : pot_fns_[5].fns)
    p_fns.push_back(fn);
  std::vector<Basis *> q_fns;
  for (Basis*& fn : pot_fns_[6].fns)
    q_fns.push_back(fn);
  std::vector<Basis *> p2_fns;
  for (Basis*& fn : pot_fns_[7].fns)
    p2_fns.push_back(fn);
  std::vector<Basis *> q2_fns;
  for (Basis*& fn : pot_fns_[8].fns)
    q2_fns.push_back(fn);

  // Set up constraint error (error from density going out of bounds of embedding function)
  Vector<double> constraint_err(mmz->config->ncells,0.0);

  // Loop over all atoms in atomvec
  for (Atom*& atom_i_ptr : mmz->atomvec->atoms) {
    // Make temporary atom and cell
    AtomJMEAM2 &atom_i = *(static_cast<AtomJMEAM2 *>(atom_i_ptr));
    Cell &cell = *mmz->config->cells[atom_i.cell_idx];

    double rho_val = 0.0; // initialize density for this atom
    double dF = 0.0;      // initialize gradient of embedding fn for this atom

    // Loop over pairs for this atom
    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairJMEAM2 &pair_ij = *(static_cast<PairJMEAM2 *>(pair_ij_ptr));  // tmp pair

      // Check that neighbor length lies in pair potential radius
      if (pair_ij.phi_knot != -1) {
        AtomJMEAM2 &atom_j = *(static_cast<AtomJMEAM2 *>(pair_ij.neigh));  // tmp atom

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

      // Check that neighbor length lies in p- potential radius
      if (pair_ij.p_knot != -1) {
        pair_ij.p = p_fns[pair_ij.p_idx]->eval_comb(pair_ij.r, &pair_ij.dp);
      } else {
        pair_ij.p = 0.0;
        pair_ij.dp = 0.0;
      } // END IF STMNT: PAIR LIES INSIDE CUTOFF FOR p2- POTENTIAL

      // Check that neighbor length lies in p2- potential radius
      if (pair_ij.p2_knot != -1) {
        pair_ij.p2 = p2_fns[pair_ij.p2_idx]->eval_comb(pair_ij.r, &pair_ij.dp2);
      } else {
        pair_ij.p2 = 0.0;
        pair_ij.dp2 = 0.0;
      } // END IF STMNT: PAIR LIES INSIDE CUTOFF FOR p2- POTENTIAL
    } // END LOOP OVER PAIRS

    // Loop over every angle formed by pairs called triplets
    for (Triplet*& triplet_ijk_ptr : atom_i.triplets) {
      TripletJMEAM2 &triplet_ijk = *(static_cast<TripletJMEAM2 *>(triplet_ijk_ptr));  // tmp triplet
      PairJMEAM2 &pair_ij = *(static_cast<PairJMEAM2 *>(triplet_ijk.pair_ij));  // tmp pairs
      PairJMEAM2 &pair_ik = *(static_cast<PairJMEAM2 *>(triplet_ijk.pair_ik));

      // The cos(theta) should always lie inside -1 ... 1
      // So store the g and g' without checking bounds
      triplet_ijk.g = g_fns[triplet_ijk.g_idx]->eval_comb(triplet_ijk.cos, &triplet_ijk.dg);
      triplet_ijk.q = q_fns[triplet_ijk.q_idx]->eval_comb(triplet_ijk.cos, &triplet_ijk.dq);
      triplet_ijk.q2 = q2_fns[triplet_ijk.q2_idx]->eval_comb(triplet_ijk.cos, &triplet_ijk.dq2);

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
    cell.energy += F_fns[atom_i.F_idx]->eval_comb(rho_val, &dF);

    // Loop over pairs for this atom to compute EAM force
    for (Pair*& pair_ij_ptr : atom_i.pairs) {
      PairJMEAM2 &pair_ij = *(static_cast<PairJMEAM2 *>(pair_ij_ptr));  // tmp pair
      AtomJMEAM2 &atom_j = *(static_cast<AtomJMEAM2 *>(pair_ij.neigh));  // tmp atom

      Vect tmp_force = pair_ij.dist * pair_ij.drho * dF;  // compute tmp force values
      atom_i.force += tmp_force;  // add in force on atom i from atom j
      atom_j.force -= tmp_force;  // subtract off force on atom j from atom i (Newton's law: action = -reaction)

      // Compute stress on cell
      tmp_force *= pair_ij.r;
      cell.stress -= pair_ij.dist & tmp_force;
    } // END 2nd LOOP OVER PAIRS

    // Loop over every angle formed by pairs called triplets
    for (Triplet*& triplet_ijk_ptr : atom_i.triplets) {
      TripletJMEAM2 &triplet_ijk = *(static_cast<TripletJMEAM2 *>(triplet_ijk_ptr));  // tmp triplet
      PairJMEAM2 &pair_ij = *(static_cast<PairJMEAM2 *>(triplet_ijk.pair_ij));  // tmp pairs
      PairJMEAM2 &pair_ik = *(static_cast<PairJMEAM2 *>(triplet_ijk.pair_ik));
      AtomJMEAM2 &atom_j = *(static_cast<AtomJMEAM2 *>(pair_ij.neigh));  // tmp atoms
      AtomJMEAM2 &atom_k = *(static_cast<AtomJMEAM2 *>(pair_ik.neigh));

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
   write punishment/constraint errors (none for this pot type)
------------------------------------------------------------------------- */

void PotJMEAM2::write_punish(const ErrorVec& e_vec, String filename, int num) const
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

  punish_fs << std::endl << "Dummy constraints (|p_max(r)| = 1)" << std::endl;
  punish_fs << std::setw(15) << std::right << "punish error"
            << std::setw(15) << std::right << "|p_max(r)|"
            << std::endl;

  double max_p_mag = std::abs(pot_fns_[5].get_max_y_mag());

  punish_fs << std::fixed
            << std::setw(15) << std::right << std::setprecision(6) << e_vec[punish_idx]*e_vec[punish_idx]
            << std::setw(15) << std::right << std::setprecision(6) << std::abs(max_p_mag)
            << std::endl;
  ++punish_idx;

  punish_fs << std::endl << "Dummy constraints (|p2_max(r)| = 1)" << std::endl;
  punish_fs << std::setw(15) << std::right << "punish error"
            << std::setw(15) << std::right << "|p2_max(r)|"
            << std::endl;

  double max_p2_mag = std::abs(pot_fns_[7].get_max_y_mag());

  punish_fs << std::fixed
            << std::setw(15) << std::right << std::setprecision(6) << e_vec[punish_idx]*e_vec[punish_idx]
            << std::setw(15) << std::right << std::setprecision(6) << std::abs(max_p2_mag)
            << std::endl;
  ++punish_idx;

  punish_fs.close();

  return;
}

/* ----------------------------------------------------------------------
  rescale MEAM part of the potential: f*f*g
------------------------------------------------------------------------- */

int PotJMEAM2::rescale_3body(const Comm& comm, std::ostream *out, int flag)
{
  if (!flag) return 0; // Don't run rescaling if flag is 0

  PotJMEAM::rescale_3body(comm, out, flag); // rescale MEAM piece first (always returns 1, unless flag=0)

  // Get scaling factor
  double max_p_mag = pot_fns_[7].get_max_y_mag();
  double b = 1.0/max_p_mag;

  // Scale beta first: f*f*g = (b*f)(b*f)(g/b^2)
  // Scale f-pot first
  std::vector<Basis *> p_fns = pot_fns_[7].fns;
  for (Basis *p_fn : p_fns)
    *p_fn *= b;

  // We scale g by 1/b^2
  std::vector<Basis *> q_fns = pot_fns_[8].fns;
  for (Basis *q_fn : q_fns)
    *q_fn /= b*b;

  // Output scaling factor to screen
  if (comm.is_root() && out)
    *out << "JMEAM2 potential scaling factor " << std::fixed << b << std::endl;

  return 1;
}
