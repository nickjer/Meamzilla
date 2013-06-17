#include "pair_jmeam2.h"
#include "../../atom.h"
#include "../../pot_fns.h"
#include "../../error.h"

#include <sstream>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

PairJMEAM2::PairJMEAM2() : PairJMEAM(), p2_idx(0), p2_knot(-1), p2(0.0), dp2(0.0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

PairJMEAM2::~PairJMEAM2()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone pair
------------------------------------------------------------------------- */

PairJMEAM2* PairJMEAM2::clone() const
{
  return new PairJMEAM2(*this);
}

/* ----------------------------------------------------------------------
   Check whether the pair i,j atoms exists in potential
   inner cutoff
------------------------------------------------------------------------- */

bool PairJMEAM2::check_pair(Atom *atom, Potential *pot)
{
  int typ1 = atom->typ;
  int typ2 = neigh->typ;

  // Setup potential
  PotFns& phi = pot->at(0);
  phi_idx = phi.get_alloy_idx(typ1, typ2);
  Basis& phi_fn = *phi.fns[phi_idx];

  // Check if pair lies inside radial cutoff
  if (r < phi_fn.get_min_rcut()) {
    std::ostringstream oss;
    oss << "ERROR: Pair distance too short for potential (r=" << r << ", pot=phi[" << phi_idx << "]) - ";
    throw Error(oss.str());
  }
  if (r <= phi_fn.get_max_rcut()) phi_knot = 1;

  // Setup potential
  PotFns& rho = pot->at(1);
  rho_idx = rho.get_alloy_idx(typ1, typ2);
  Basis& rho_fn = *rho.fns[rho_idx];

  // Check if pair lies inside radial cutoff
  if (r < rho_fn.get_min_rcut()) {
    std::ostringstream oss;
    oss << "ERROR: Pair distance too short for potential (r=" << r << ", pot=rho[" << rho_idx << "]) - ";
    throw Error(oss.str());
  }
  if (r <= rho_fn.get_max_rcut()) rho_knot = 1;

  // Setup potential
  PotFns& f = pot->at(3);
  f_idx = f.get_alloy_idx(typ1, typ2);
  Basis& f_fn = *f.fns[f_idx];

  // Check if pair lies inside radial cutoff
  if (r < f_fn.get_min_rcut()) {
    std::ostringstream oss;
    oss << "ERROR: Pair distance too short for potential (r=" << r << ", pot=f[" << f_idx << "]) - ";
    throw Error(oss.str());
  }
  if (r <= f_fn.get_max_rcut()) f_knot = 1;

  // Setup potential
  PotFns& p = pot->at(5);
  p_idx = p.get_alloy_idx(typ1, typ2);
  Basis& p_fn = *p.fns[p_idx];

  // Check if pair lies inside radial cutoff
  if (r < p_fn.get_min_rcut()) {
    std::ostringstream oss;
    oss << "ERROR: Pair distance too short for potential (r=" << r << ", pot=p[" << p_idx << "]) - ";
    throw Error(oss.str());
  }
  if (r <= p_fn.get_max_rcut()) p_knot = 1;

  // Setup potential
  PotFns& p2 = pot->at(7);
  p2_idx = p2.get_alloy_idx(typ1, typ2);
  Basis& p2_fn = *p2.fns[p2_idx];

  // Check if pair lies inside radial cutoff
  if (r < p2_fn.get_min_rcut()) {
    std::ostringstream oss;
    oss << "ERROR: Pair distance too short for potential (r=" << r << ", pot=p2[" << p2_idx << "]) - ";
    throw Error(oss.str());
  }
  if (r <= p2_fn.get_max_rcut()) p2_knot = 1;

  // If pair lies outside of all potentials then return false
  if (phi_knot==-1 && rho_knot==-1 && f_knot==-1 && p_knot==-1 && p2_knot==-1)
    return false;

  return true;
}
