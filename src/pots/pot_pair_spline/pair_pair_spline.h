/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(pair/spline,PairPairSpline)

#else

#ifndef MEAMZ_PAIR_PAIR_SPLINE_H
#define MEAMZ_PAIR_PAIR_SPLINE_H

#include "../pot_pair/pair_pair.h"

namespace MEAMZ_NS
{

class PairPairSpline : public PairPair
{
public:
  double phi_shift;

  PairPairSpline();
  virtual ~PairPairSpline();

  virtual PairPairSpline* clone() const;          // Copy polymorphic objects

  virtual bool check_pair(Atom *, Potential *);   // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_PAIR_PAIR_SPLINE_H
#endif // PAIR_CLASS
