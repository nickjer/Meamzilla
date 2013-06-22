/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "pair.h"


using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Pair::Pair() : neigh(nullptr), r(0.0), invr(0.0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

Pair::~Pair()
{
  //dtor
}

/* ----------------------------------------------------------------------
   setup pair with user specified conditions
------------------------------------------------------------------------- */

void Pair::setup_pair(Atom *_neigh, Vect _total_dist)
{
  neigh = _neigh;
  r = mag(_total_dist);
  invr = 1.0/r;
  dist = _total_dist/r;
  return;
}
