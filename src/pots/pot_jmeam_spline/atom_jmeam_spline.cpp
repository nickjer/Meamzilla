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
