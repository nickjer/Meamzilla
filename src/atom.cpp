/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "atom.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Atom::Atom() : typ(0), cell_idx(0), atom_idx(0), global_idx(0), npairs(0), ntriplets(0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

Atom::~Atom()
{
  for (Pair*& pair_ptr : pairs)
    if (pair_ptr) delete pair_ptr;
  for (Triplet*& triplet_ptr : triplets)
    if (triplet_ptr) delete triplet_ptr;
}

/* ----------------------------------------------------------------------
   setup atom with user specified conditions
------------------------------------------------------------------------- */

void Atom::setup_atom(int _idx, int _typ, Vect _pos, Vect _force)
{
  typ = _typ;
  pos = _pos;
  force0 = _force;
  atom_idx = _idx;

  setup_atom_pot(); // setup potential specific atom properties
  return;
}

/* ----------------------------------------------------------------------
   setup potential specific atom properties
------------------------------------------------------------------------- */

void Atom::setup_atom_pot()
{
  return;
}


/* ----------------------------------------------------------------------
   add this pair to list of pairs
------------------------------------------------------------------------- */

void Atom::push_back(Pair *pair_ptr)
{
  pairs.push_back(pair_ptr);
  ++npairs;
}

/* ----------------------------------------------------------------------
   add this triplet to list of triplets
------------------------------------------------------------------------- */

void Atom::push_back(Triplet *triplet_ptr)
{
  triplets.push_back(triplet_ptr);
  ++ntriplets;
}
