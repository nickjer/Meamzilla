/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(jmeam2,PairJMEAM2)

#else

#ifndef MEAMZ_PAIR_JMEAM2_H
#define MEAMZ_PAIR_JMEAM2_H

#include "../pot_jmeam/pair_jmeam.h"

namespace MEAMZ_NS
{

class PairJMEAM2 : public PairJMEAM
{
public:
  int p2_idx;                        // Index for p- alloy potential used for this pair
  int p2_knot;                       // -1 if outside range of p- potential, anything else otherwise

  double p2;                         // p-potential value for pair
  double dp2;                        // Derivative of p- potential for this pair: d(p)/dr

  PairJMEAM2();
  virtual ~PairJMEAM2();

  virtual PairJMEAM2* clone() const;            // Copy polymorphic objects

  virtual bool check_pair(Atom *, Potential *); // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_PAIR_JMEAM2_H
#endif // PAIR_CLASS
