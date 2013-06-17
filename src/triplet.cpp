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
