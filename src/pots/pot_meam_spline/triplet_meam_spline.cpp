/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "triplet_meam_spline.h"
#include "../../atom.h"
#include "../../pot_fns.h"
#include "../../spline.h"

#include <sstream>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

TripletMEAMSpline::TripletMEAMSpline() : TripletMEAM()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

TripletMEAMSpline::~TripletMEAMSpline()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone pair
------------------------------------------------------------------------- */

TripletMEAMSpline* TripletMEAMSpline::clone() const
{
  return new TripletMEAMSpline(*this);
}

/* ----------------------------------------------------------------------
   Check whether the pair i,j atoms exists in potential
   inner cutoff
------------------------------------------------------------------------- */

bool TripletMEAMSpline::check_triplet(Atom *atom, Potential *pot)
{
  int typ_i = atom->typ;
  int typ_j = pair_ij->neigh->typ;
  int typ_k = pair_ik->neigh->typ;

  // Setup potentials
  PotFns& f = pot->at(3);
  Basis& f_ij_fn = *f.fns[ f.get_2body_alloy_idx(typ_i, typ_j) ];
  Basis& f_ik_fn = *f.fns[ f.get_2body_alloy_idx(typ_i, typ_k) ];

  PotFns& g = pot->at(4);
  g_idx = g.get_3body_alloy_idx(typ_i, typ_j, typ_k);  // only depends on origin atom_i's type
  Spline& g_spline = *static_cast<Spline*>(g.fns[g_idx]);

  // Check if both pairs lie inside radial cutoffs of respective fn's
  if (pair_ij->r > f_ij_fn.get_max_rcut() || pair_ik->r > f_ik_fn.get_max_rcut()) g_knot = -1;
  g_knot = g_spline.get_knotshift( cos, g_shift );

  // If pair lies outside of all potentials then return false
  if (g_knot==-1)
    return false;

  return true;
}
