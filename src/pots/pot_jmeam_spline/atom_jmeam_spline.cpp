/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "atom_jmeam_spline.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

AtomJMEAMSpline::AtomJMEAMSpline() : AtomJMEAM()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

AtomJMEAMSpline::~AtomJMEAMSpline()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone atom
------------------------------------------------------------------------- */

AtomJMEAMSpline* AtomJMEAMSpline::clone() const
{
  return new AtomJMEAMSpline(*this);
}
