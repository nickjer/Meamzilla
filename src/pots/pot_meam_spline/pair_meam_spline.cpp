/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "pair_meam_spline.h"
#include "../../atom.h"
#include "../../pot_fns.h"
#include "../../spline.h"
#include "../../error.h"

#include <sstream>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

PairMEAMSpline::PairMEAMSpline() : PairMEAM()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

PairMEAMSpline::~PairMEAMSpline()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone pair
------------------------------------------------------------------------- */

PairMEAMSpline* PairMEAMSpline::clone() const
{
  return new PairMEAMSpline(*this);
}

/* ----------------------------------------------------------------------
   Check whether the pair i,j atoms exists in potential
   inner cutoff
------------------------------------------------------------------------- */

bool PairMEAMSpline::check_pair(Atom *atom, Potential *pot)
{
  int typ1 = atom->typ;
  int typ2 = neigh->typ;

  // Setup potential
  PotFns& phi = pot->at(0);
  phi_idx = phi.get_2body_alloy_idx(typ1, typ2);
  Spline& phi_spline = *static_cast<Spline*>(phi.fns[phi_idx]);

  // Check if pair lies inside radial cutoff
  if (r < phi_spline.get_min_rcut()) {
    std::ostringstream oss;
    oss << "ERROR: Pair distance too short for potential (r=" << r << ", pot=phi[" << phi_idx << "]) - ";
    throw Error(oss.str());
  }
  phi_knot = phi_spline.get_knotshift( r, phi_shift );

  // Setup potential
  PotFns& rho = pot->at(1);
  rho_idx = rho.get_2body_alloy_idx(typ1, typ2);
  Spline& rho_spline = *static_cast<Spline*>(rho.fns[rho_idx]);

  // Check if pair lies inside radial cutoff
  if (r < rho_spline.get_min_rcut()) {
    std::ostringstream oss;
    oss << "ERROR: Pair distance too short for potential (r=" << r << ", pot=rho[" << rho_idx << "]) - ";
    throw Error(oss.str());
  }
  rho_knot = rho_spline.get_knotshift( r, rho_shift );

  // Setup potential
  PotFns& f = pot->at(3);
  f_idx = f.get_2body_alloy_idx(typ1, typ2);
  Spline& f_spline = *static_cast<Spline*>(f.fns[f_idx]);

  // Check if pair lies inside radial cutoff
  if (r < f_spline.get_min_rcut()) {
    std::ostringstream oss;
    oss << "ERROR: Pair distance too short for potential (r=" << r << ", pot=f[" << f_idx << "]) - ";
    throw Error(oss.str());
  }
  f_knot = f_spline.get_knotshift( r, f_shift );

  // If pair lies outside of all potentials then return false
  if (phi_knot==-1 && rho_knot==-1 && f_knot==-1)
    return false;

  return true;
}
