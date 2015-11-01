/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifndef MEAMZ_POT_FNS_H
#define MEAMZ_POT_FNS_H

#include "basis.h"

namespace MEAMZ_NS
{

class PotFns
{
public:
  enum class AlloyType;               // Class defining what this potential depends on in alloy form
  enum class FnType;                  // Class defining what type of fn it is (e.g., radial)

  int n_fns;                          // Number of alloy potentials (read below for # based on alloy type)
  std::vector<Basis *> fns;           // Each potential type has a list of alloy dependent potentials
                                      // (e.g., phi_ij = phi_AA, phi_AB, phi_BB for 2 atom type alloy)

  PotFns();
  PotFns(const PotFns&);
  virtual ~PotFns();

  PotFns& operator=(const PotFns&);

  int is_radial() const;                    // Returns whether this potential function is radial or not
  double get_max_rcut() const;              // Get maximum radial cutoff from all potentials
  double get_max_y_mag() const;             // Get maximum |y| value from all fns in this potential
  int get_1body_alloy_idx(int) const;       // Get index of 1-body alloy basis fn knowing types of atoms
  int get_2body_alloy_idx(int, int) const;       // Get index of 2-body alloy basis fn knowing types of atoms
  int get_3body_alloy_idx(int, int, int) const;  // Get index of 3-body alloy basis fn knowing types of atoms
  String get_common_basis_type() const;     // Get the basis type common among all basis fns for this pot fn (return "undefined" if no common type)
  int get_ncoeff() const;                   // Get number of adjustable coefficients in basis fns
  double& coeff(int&);                      // Return value of specified coeffient
  void modify_coeff(int, double, double);   // Perform modification to fn at specified coefficient

  void check_if_similar(PotFns&) const;     // Check that a pot fn is similar enough to this
                                            // pot fn for pairs/triplets to have identical info
  void refresh_basis();                     // Refresh values in basis set to account for changes made earlier


  // Initialize alloy fns with correct fn type, alloy type, and # of fns
  void setup_pot_fns(int, AlloyType, FnType);

  int read_pot_fns(const StringList&, int); // Read in each alloy potential for this unique potential
  void write_pot_fns(std::ostream&) const;  // Write out potential fn for each alloy fn
  void write_lmp_fns(std::ostream&) const;  // Write out potential fn for each alloy fn in lammps format

  void communicate(const Comm&, int);       // Communicate data to all procs within group

protected:

private:
  AlloyType alloy_type_;         // Each potential fn will depend on alloy types
                                 // See:
                                 //     PotFns::AlloyType

  FnType fn_type_;               // Whether this potential fn is a radial fn or not

  int ntypes_;                   // Number of atom types
};

// Define the alloy dependency for this alloy potential
enum class PotFns::AlloyType
{
  Atom_i,       // depends on atom in question, atom_i (A, B)
  Pair_i,       // depends on atom in question for pair fn, atom_i (A, B)
  Pair_j,       // depends on neighbor atom for pair fn, atom_j (A, B)
  Pair_ij,      // depends on pair of atoms i and j (AA, AB, BB)
  Triplet_i,    // depends on atom in question for triplet fn, atom_i (A, B)
  Triplet_j,    // depends on neighbor atom for triplet fn, atom_j (A, B)
  Triplet_k,    // depends on neighbor atom for triplet fn, atom_k (A, B)
  Triplet_i_jk, // depends on atom and combination of pairs (AAA, AAB, ABB, BAA, BAB, BBB)
  Undefined     // unknown and will kill the program
};

// Define the functional form for this alloy potential
enum class PotFns::FnType
{
  Radial,     // radial fn
  Undefined   // undefined fn
};

} // namespace MEAMZ_NS

#endif // MEAMZ_POT_FNS_H
