/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "atom_meam.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

AtomMEAM::AtomMEAM() : AtomEAM()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

AtomMEAM::~AtomMEAM()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone atom
------------------------------------------------------------------------- */

AtomMEAM* AtomMEAM::clone() const
{
  return new AtomMEAM(*this);
}

/* ----------------------------------------------------------------------
   value used when sorting atoms
------------------------------------------------------------------------- */

int AtomMEAM::cmp_value() const
{
  return ntriplets;
}
