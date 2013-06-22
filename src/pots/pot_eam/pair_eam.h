/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(eam,PairEAM)

#else

#ifndef MEAMZ_PAIR_EAM_H
#define MEAMZ_PAIR_EAM_H

#include "../pot_pair/pair_pair.h"

namespace MEAMZ_NS
{

class PairEAM : public PairPair
{
public:
                                    // Carry over phi info from PairPair
  int rho_idx;                      // Index for rho alloy potential used for this pair
  int rho_knot;                     // -1 if outside range of rho potential, anything else otherwise

  double drho;                      // Derivative of rho potential for this pair: d(rho)/dr

  PairEAM();
  virtual ~PairEAM();

  virtual PairEAM* clone() const;               // Copy polymorphic objects

  virtual bool check_pair(Atom *, Potential *); // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_PAIR_EAM_H
#endif // PAIR_CLASS
