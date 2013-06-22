/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "triplet_jmeam2.h"
#include "../../atom.h"
#include "../../pot_fns.h"

#include <sstream>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

TripletJMEAM2::TripletJMEAM2() : TripletJMEAM(), q2_idx(0), q2_knot(-1), q2(0.0), dq2(0.0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

TripletJMEAM2::~TripletJMEAM2()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone pair
------------------------------------------------------------------------- */

TripletJMEAM2* TripletJMEAM2::clone() const
{
  return new TripletJMEAM2(*this);
}

/* ----------------------------------------------------------------------
   Check whether the pair i,j atoms exists in potential
   inner cutoff
------------------------------------------------------------------------- */

bool TripletJMEAM2::check_triplet(Atom *atom, Potential *pot)
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
  //Basis& g_fn = *g.fns[g_idx];

  // Check if both pairs lie inside radial cutoffs of respective fn's
  if (pair_ij->r > f_ij_fn.get_max_rcut() || pair_ik->r > f_ik_fn.get_max_rcut()) g_knot = -1;
  else g_knot = 1;

  // Setup potentials
  PotFns& p = pot->at(5);
  Basis& p_ij_fn = *p.fns[ p.get_alloy_idx(typ_i, typ_j) ];
  Basis& p_ik_fn = *p.fns[ p.get_alloy_idx(typ_i, typ_k) ];

  PotFns& q = pot->at(6);
  q_idx = q.get_alloy_idx(typ_i, 0);  // only depends on origin atom_i's type
  //Basis& q_fn = *q.fns[q_idx];

  // Check if both pairs lie inside radial cutoffs of respective fn's
  if (pair_ij->r > p_ij_fn.get_max_rcut() || pair_ik->r > p_ik_fn.get_max_rcut()) q_knot = -1;
  else q_knot = 1;

  // Setup potentials
  PotFns& p2 = pot->at(7);
  Basis& p2_ij_fn = *p2.fns[ p2.get_alloy_idx(typ_i, typ_j) ];
  Basis& p2_ik_fn = *p2.fns[ p2.get_alloy_idx(typ_i, typ_k) ];

  PotFns& q2 = pot->at(8);
  q2_idx = q2.get_alloy_idx(typ_i, 0);  // only depends on origin atom_i's type
  //Basis& q2_fn = *q2.fns[q_idx];

  // Check if both pairs lie inside radial cutoffs of respective fn's
  if (pair_ij->r > p2_ij_fn.get_max_rcut() || pair_ik->r > p2_ik_fn.get_max_rcut()) q2_knot = -1;
  else q2_knot = 1;

  // If pair lies outside of all potentials then return false
  if (g_knot==-1 && q_knot==-1 && q2_knot==-1)
    return false;

  return true;
}
