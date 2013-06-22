/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "triplet.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Triplet::Triplet() : pair_ij(nullptr), pair_ik(nullptr), cos(0.0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

Triplet::~Triplet()
{
  //dtor
}

/* ----------------------------------------------------------------------
   setup pair with user specified conditions
------------------------------------------------------------------------- */

void Triplet::setup_triplet(Pair *_pij, Pair *_pik)
{
  pair_ij = _pij;
  pair_ik = _pik;
  cos = _pij->dist * _pik->dist;
  return;
}
