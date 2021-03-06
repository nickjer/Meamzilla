/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(pair,PairPair)

#else

#ifndef MEAMZ_PAIR_PAIR_H
#define MEAMZ_PAIR_PAIR_H

#include "../../pair.h"

namespace MEAMZ_NS
{

class PairPair : public Pair
{
public:
  int phi_idx;                      // Index for phi alloy potential used for this pair
  int phi_knot;                     // -1 if outside range of phi potential, anything else otherwise

  PairPair();
  virtual ~PairPair();

  virtual PairPair* clone() const;              // Copy polymorphic objects

  virtual bool check_pair(Atom *, Potential *); // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_PAIR_PAIR_H
#endif // PAIR_CLASS
