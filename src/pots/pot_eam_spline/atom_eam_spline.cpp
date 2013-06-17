#include "atom_eam_spline.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

AtomEAMSpline::AtomEAMSpline() : AtomEAM()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

AtomEAMSpline::~AtomEAMSpline()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone atom
------------------------------------------------------------------------- */

AtomEAMSpline* AtomEAMSpline::clone() const
{
  return new AtomEAMSpline(*this);
}
