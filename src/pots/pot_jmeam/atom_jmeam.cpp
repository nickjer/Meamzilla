/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

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
