#include "atom_pair.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

AtomPair::AtomPair() : Atom()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

AtomPair::~AtomPair()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone atom
------------------------------------------------------------------------- */

AtomPair* AtomPair::clone() const
{
  return new AtomPair(*this);
}

/* ----------------------------------------------------------------------
   value used when sorting atoms
------------------------------------------------------------------------- */

int AtomPair::cmp_value() const
{
  return npairs;
}
