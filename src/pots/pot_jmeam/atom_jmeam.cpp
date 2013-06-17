#include "atom_jmeam.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

AtomJMEAM::AtomJMEAM() : AtomMEAM()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

AtomJMEAM::~AtomJMEAM()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone atom
------------------------------------------------------------------------- */

AtomJMEAM* AtomJMEAM::clone() const
{
  return new AtomJMEAM(*this);
}
