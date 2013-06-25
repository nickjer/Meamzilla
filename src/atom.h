/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifndef MEAMZ_ATOM_H
#define MEAMZ_ATOM_H

#include "vect.h"
#include "pair.h"
#include "triplet.h"
#include <vector>

namespace MEAMZ_NS
{

class Atom
{
public:
  int typ;       // Atom type
  Vect pos;      // 3D position vector
  Vect force;    // 3D computed force vector
  Vect force0;   // 3D real force vector

  int cell_idx;      // Global index of supercell this atom is in
  int atom_idx;      // Index of atom in supercell
  int global_idx;    // Global atom index

  // Pairs for this atom
  int npairs;
  std::vector<Pair *> pairs;

  int ntriplets;
  std::vector<Triplet *> triplets;

  Atom();
  virtual ~Atom();

  virtual Atom* clone() const = 0;        // Copy polymorphic objects

  virtual int cmp_value() const = 0;      // Value used when sorting atoms
  void setup_atom(int, int, Vect, Vect);  // Setup atom with user specified conditions
  virtual void setup_atom_pot();          // Setup potential specific atom properties
  void push_back(Pair *);                 // Add this pair to list of pairs
  void push_back(Triplet *);              // Add this triplet to list of triplets

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_H
