/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(jmeam2/spline,AtomJMEAM2Spline)

#else

#ifndef MEAMZ_ATOM_JMEAM2_SPLINE_H
#define MEAMZ_ATOM_JMEAM2_SPLINE_H

#include "../pot_jmeam2/atom_jmeam2.h"

namespace MEAMZ_NS
{

class AtomJMEAM2Spline : public AtomJMEAM2
{
public:
  AtomJMEAM2Spline();
  virtual ~AtomJMEAM2Spline();

  virtual AtomJMEAM2Spline* clone() const;  // Copy polymorphic objects

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_JMEAM2_SPLINE_H
#endif // ATOM_CLASS
