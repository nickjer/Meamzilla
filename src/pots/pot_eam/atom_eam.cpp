/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "atom_eam.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

AtomEAM::AtomEAM() : AtomPair(), F_idx(0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

AtomEAM::~AtomEAM()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone atom
------------------------------------------------------------------------- */

AtomEAM* AtomEAM::clone() const
{
  return new AtomEAM(*this);
}

/* ----------------------------------------------------------------------
   setup potential specific atom properties
------------------------------------------------------------------------- */

void AtomEAM::setup_atom_pot()
{
  F_idx = typ;  // embedding fn depends on origin atom_i's type
  return;
}

