#include "atom_pair_spline.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

AtomPairSpline::AtomPairSpline() : AtomPair()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

AtomPairSpline::~AtomPairSpline()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone atom
------------------------------------------------------------------------- */

AtomPairSpline* AtomPairSpline::clone() const
{
  return new AtomPairSpline(*this);
}
