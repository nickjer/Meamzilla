/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(jmeam,PairJMEAM)

#else

#ifndef MEAMZ_PAIR_JMEAM_H
#define MEAMZ_PAIR_JMEAM_H

#include "../pot_meam/pair_meam.h"

namespace MEAMZ_NS
{

class PairJMEAM : public PairMEAM
{
public:
  int p_idx;                        // Index for p- alloy potential used for this pair
  int p_knot;                       // -1 if outside range of p- potential, anything else otherwise

  double p;                         // p-potential value for pair
  double dp;                        // Derivative of p- potential for this pair: d(p)/dr

  PairJMEAM();
  virtual ~PairJMEAM();

  virtual PairJMEAM* clone() const;             // Copy polymorphic objects

  virtual bool check_pair(Atom *, Potential *); // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_PAIR_JMEAM_H
#endif // PAIR_CLASS
