/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "triplet_jmeam_spline.h"
#include "../../atom.h"
#include "../../pot_fns.h"
#include "../../spline.h"

#include <sstream>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

TripletJMEAMSpline::TripletJMEAMSpline() : TripletJMEAM()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

TripletJMEAMSpline::~TripletJMEAMSpline()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone pair
------------------------------------------------------------------------- */

TripletJMEAMSpline* TripletJMEAMSpline::clone() const
{
  return new TripletJMEAMSpline(*this);
}

/* ----------------------------------------------------------------------
   Check whether the pair i,j atoms exists in potential
   inner cutoff
------------------------------------------------------------------------- */

bool TripletJMEAMSpline::check_triplet(Atom *atom, Potential *pot)
{
  int typ_i = atom->typ;
  int typ_j = pair_ij->neigh->typ;
  int typ_k = pair_ik->neigh->typ;

  // Setup potentials
  PotFns& f = pot->at(3);
  Basis& f_ij_fn = *f.fns[ f.get_alloy_idx(typ_i, typ_j) ];
  Basis& f_ik_fn = *f.fns[ f.get_alloy_idx(typ_i, typ_k) ];

  PotFns& g = pot->at(4);
  g_idx = g.get_alloy_idx(typ_i, 0);  // only depends on origin atom_i's type
  Spline& g_spline = *static_cast<Spline*>(g.fns[g_idx]);

  // Check if both pairs lie inside radial cutoffs of respective fn's
  if (pair_ij->r > f_ij_fn.get_max_rcut() || pair_ik->r > f_ik_fn.get_max_rcut()) g_knot = -1;
  g_knot = g_spline.get_knotshift( cos, g_shift );

  // Setup potentials
  PotFns& p = pot->at(5);
  Basis& p_ij_fn = *p.fns[ p.get_alloy_idx(typ_i, typ_j) ];
  Basis& p_ik_fn = *p.fns[ p.get_alloy_idx(typ_i, typ_k) ];

  PotFns& q = pot->at(6);
  q_idx = q.get_alloy_idx(typ_i, 0);  // only depends on origin atom_i's type
  Spline& q_spline = *static_cast<Spline*>(q.fns[q_idx]);

  // Check if both pairs lie inside radial cutoffs of respective fn's
  if (pair_ij->r > p_ij_fn.get_max_rcut() || pair_ik->r > p_ik_fn.get_max_rcut()) q_knot = -1;
  q_knot = q_spline.get_knotshift( cos, q_shift );

  // If pair lies outside of all potentials then return false
  if (g_knot==-1 && q_knot==-1)
    return false;

  return true;
}
