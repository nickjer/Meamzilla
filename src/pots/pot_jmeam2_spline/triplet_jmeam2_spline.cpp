/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "triplet_jmeam2_spline.h"
#include "../../atom.h"
#include "../../pot_fns.h"
#include "../../spline.h"

#include <sstream>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

TripletJMEAM2Spline::TripletJMEAM2Spline() : TripletJMEAM2()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

TripletJMEAM2Spline::~TripletJMEAM2Spline()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone pair
------------------------------------------------------------------------- */

TripletJMEAM2Spline* TripletJMEAM2Spline::clone() const
{
  return new TripletJMEAM2Spline(*this);
}

/* ----------------------------------------------------------------------
   Check whether the pair i,j atoms exists in potential
   inner cutoff
------------------------------------------------------------------------- */

bool TripletJMEAM2Spline::check_triplet(Atom *atom, Potential *pot)
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

  // Setup potentials
  PotFns& p = pot->at(5);
  Basis& p_ij_fn = *p.fns[ p.get_2body_alloy_idx(typ_i, typ_j) ];
  Basis& p_ik_fn = *p.fns[ p.get_2body_alloy_idx(typ_i, typ_k) ];

  PotFns& q = pot->at(6);
  q_idx = q.get_3body_alloy_idx(typ_i, typ_j, typ_k);  // only depends on origin atom_i's type
  Spline& q_spline = *static_cast<Spline*>(q.fns[q_idx]);

  // Check if both pairs lie inside radial cutoffs of respective fn's
  if (pair_ij->r > p_ij_fn.get_max_rcut() || pair_ik->r > p_ik_fn.get_max_rcut()) q_knot = -1;
  q_knot = q_spline.get_knotshift( cos, q_shift );

  // Setup potentials
  PotFns& p2 = pot->at(7);
  Basis& p2_ij_fn = *p2.fns[ p2.get_2body_alloy_idx(typ_i, typ_j) ];
  Basis& p2_ik_fn = *p2.fns[ p2.get_2body_alloy_idx(typ_i, typ_k) ];

  PotFns& q2 = pot->at(8);
  q2_idx = q2.get_3body_alloy_idx(typ_i, typ_j, typ_k);  // only depends on origin atom_i's type
  Spline& q2_spline = *static_cast<Spline*>(q2.fns[q2_idx]);

  // Check if both pairs lie inside radial cutoffs of respective fn's
  if (pair_ij->r > p2_ij_fn.get_max_rcut() || pair_ik->r > p2_ik_fn.get_max_rcut()) q2_knot = -1;
  q2_knot = q2_spline.get_knotshift( cos, q2_shift );

  // If pair lies outside of all potentials then return false
  if (g_knot==-1 && q_knot==-1 && q2_knot==-1)
    return false;

  return true;
}
