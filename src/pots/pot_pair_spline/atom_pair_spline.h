/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(pair/spline,AtomPairSpline)

#else

#ifndef MEAMZ_ATOM_PAIR_SPLINE_H
#define MEAMZ_ATOM_PAIR_SPLINE_H

#include "../pot_pair/atom_pair.h"

namespace MEAMZ_NS
{

class AtomPairSpline : public AtomPair
{
public:
  AtomPairSpline();
  virtual ~AtomPairSpline();

  virtual AtomPairSpline* clone() const;  // Copy polymorphic objects

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_PAIR_SPLINE_H
#endif // ATOM_CLASS
