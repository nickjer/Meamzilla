/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(meam,PairMEAM)

#else

#ifndef MEAMZ_PAIR_MEAM_H
#define MEAMZ_PAIR_MEAM_H

#include "../pot_eam/pair_eam.h"

namespace MEAMZ_NS
{

class PairMEAM : public PairEAM
{
public:
                                    // Carry over phi/rho info from PairEAM
  int f_idx;                        // Index for f- alloy potential used for this pair
  int f_knot;                       // -1 if outside range of f- potential, anything else otherwise

  double f;                         // f-potential value for pair
  double df;                        // Derivative of f- potential for this pair: d(f)/dr

  PairMEAM();
  virtual ~PairMEAM();

  virtual PairMEAM* clone() const;              // Copy polymorphic objects

  virtual bool check_pair(Atom *, Potential *); // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_PAIR_MEAM_H
#endif // PAIR_CLASS
