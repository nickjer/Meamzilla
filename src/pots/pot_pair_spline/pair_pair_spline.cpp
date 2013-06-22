/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "pair_pair_spline.h"
#include "../../atom.h"
#include "../../pot_fns.h"
#include "../../spline.h"
#include "../../error.h"

#include <sstream>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

PairPairSpline::PairPairSpline() : PairPair()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

PairPairSpline::~PairPairSpline()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone pair
------------------------------------------------------------------------- */

PairPairSpline* PairPairSpline::clone() const
{
  return new PairPairSpline(*this);
}

/* ----------------------------------------------------------------------
   Check whether the pair i,j atoms exists in potential
   inner cutoff, and also set the appropriate
   potential indices
------------------------------------------------------------------------- */

bool PairPairSpline::check_pair(Atom *atom, Potential *pot)
{
  int typ1 = atom->typ;
  int typ2 = neigh->typ;

  // Setup spline potential
  PotFns& phi = pot->at(0);
  phi_idx = phi.get_alloy_idx(typ1, typ2);
  Spline& phi_spline = *static_cast<Spline*>(phi.fns[phi_idx]);

  // Get lower bound knot and shift for this 'r' value
  if (r < phi_spline.get_min_rcut()) {
    std::ostringstream oss;
    oss << "ERROR: Pair distance too short for potential (r=" << r << ", pot=phi[" << phi_idx << "]) - ";
    throw Error(oss.str());
  }
  phi_knot = phi_spline.get_knotshift( r, phi_shift );

  // If pair lies outside of all potentials then return false
  if (phi_knot==-1)
    return false;

  return true;
}
