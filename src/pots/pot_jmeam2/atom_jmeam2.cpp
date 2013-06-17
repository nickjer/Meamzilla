#include "atom_jmeam2.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

AtomJMEAM2::AtomJMEAM2() : AtomJMEAM()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

AtomJMEAM2::~AtomJMEAM2()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone atom
------------------------------------------------------------------------- */

AtomJMEAM2* AtomJMEAM2::clone() const
{
  return new AtomJMEAM2(*this);
}
