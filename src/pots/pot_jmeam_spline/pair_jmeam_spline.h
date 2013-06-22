/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(jmeam/spline,PairJMEAMSpline)

#else

#ifndef MEAMZ_PAIR_JMEAM_SPLINE_H
#define MEAMZ_PAIR_JMEAM_SPLINE_H

#include "../pot_jmeam/pair_jmeam.h"

namespace MEAMZ_NS
{

class PairJMEAMSpline : public PairJMEAM
{
public:
  double phi_shift;
  double rho_shift;
  double f_shift;
  double p_shift;

  PairJMEAMSpline();
  virtual ~PairJMEAMSpline();

  virtual PairJMEAMSpline* clone() const;       // Copy polymorphic objects

  virtual bool check_pair(Atom *, Potential *); // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_PAIR_JMEAM_H
#endif // PAIR_CLASS
