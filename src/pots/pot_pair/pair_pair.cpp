/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "pair_pair.h"
#include "../../atom.h"
#include "../../pot_fns.h"
#include "../../error.h"

#include <sstream>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

PairPair::PairPair() : Pair(), phi_idx(0), phi_knot(-1)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

PairPair::~PairPair()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone pair
------------------------------------------------------------------------- */

PairPair* PairPair::clone() const
{
  return new PairPair(*this);
}

/* ----------------------------------------------------------------------
   Check whether the pair i,j atoms exists in potential
   inner cutoff
------------------------------------------------------------------------- */

bool PairPair::check_pair(Atom *atom, Potential *pot)
{
  int typ1 = atom->typ;
  int typ2 = neigh->typ;

  // Setup potential
  PotFns& phi = pot->at(0);
  phi_idx = phi.get_2body_alloy_idx(typ1, typ2);
  Basis& phi_fn = *phi.fns[phi_idx];

  // Check if pair lies inside radial cutoff
  if (r < phi_fn.get_min_rcut()) {
    std::ostringstream oss;
    oss << "ERROR: Pair distance too short for potential (r=" << r << ", pot=phi[" << phi_idx << "]) - ";
    throw Error(oss.str());
  }
  if (r <= phi_fn.get_max_rcut()) phi_knot = 1;

  // If pair lies outside of all potentials then return false
  if (phi_knot==-1)
    return false;

  return true;
}
