/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifndef MEAMZ_TRIPLET_H
#define MEAMZ_TRIPLET_H

#include "pair.h"

namespace MEAMZ_NS
{

class Triplet
{
public:
  Pair *pair_ij, *pair_ik;    // Pairs that make up 3-body term
  double cos;                 // cos(theta_ijk) of 3-body term

  Triplet();
  virtual ~Triplet();

  virtual Triplet* clone() const = 0;                   // Copy polymorphic objects

  void setup_triplet(Pair *, Pair *);                   // Setup triplet with user specified conditions
  virtual bool check_triplet(Atom *, Potential *) = 0;  // Check that this triplet is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_TRIPLET_H
