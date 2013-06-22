/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "atom_jmeam2_spline.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

AtomJMEAM2Spline::AtomJMEAM2Spline() : AtomJMEAM2()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

AtomJMEAM2Spline::~AtomJMEAM2Spline()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone atom
------------------------------------------------------------------------- */

AtomJMEAM2Spline* AtomJMEAM2Spline::clone() const
{
  return new AtomJMEAM2Spline(*this);
}
