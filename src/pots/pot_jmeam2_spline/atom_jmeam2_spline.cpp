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
