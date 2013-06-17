#include "atom_meam_spline.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

AtomMEAMSpline::AtomMEAMSpline() : AtomMEAM()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

AtomMEAMSpline::~AtomMEAMSpline()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone atom
------------------------------------------------------------------------- */

AtomMEAMSpline* AtomMEAMSpline::clone() const
{
  return new AtomMEAMSpline(*this);
}
