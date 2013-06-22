/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "pot_fns.h"
#include "error.h"
#include "style_basis.h"
#include <ostream>
#include <sstream>
#include <cmath>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

PotFns::PotFns() : n_fns(0), alloy_type_(AlloyType::Undefined), fn_type_(FnType::Undefined), ntypes_(0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

PotFns::PotFns(const PotFns& pot_fn) : n_fns(pot_fn.n_fns), alloy_type_(pot_fn.alloy_type_),
                                       fn_type_(pot_fn.fn_type_), ntypes_(pot_fn.ntypes_)
{
  fns.resize(n_fns, nullptr);
  for (int i=0; i<n_fns; ++i)
    if (pot_fn.fns[i])
      fns[i] = pot_fn.fns[i]->clone();
    else
      fns[i] = nullptr;
}

/* ---------------------------------------------------------------------- */

PotFns::~PotFns()
{
  // Clean up memory
  for (Basis*& fn : fns)
    if (fn) delete fn;
}

/* ---------------------------------------------------------------------- */

PotFns& PotFns::operator=(const PotFns& pot_fn)
{
  n_fns = pot_fn.n_fns;
  alloy_type_ = pot_fn.alloy_type_;
  fn_type_ = pot_fn.fn_type_;
  ntypes_ = pot_fn.ntypes_;

  fns.resize(n_fns, nullptr);
  for (int i=0; i<n_fns; ++i) {
    if (fns[i]) delete fns[i];  // remove old fn

    // clone fn if it exists
    if (pot_fn.fns[i])
      fns[i] = pot_fn.fns[i]->clone();
    else
      fns[i] = nullptr;
  }

  return *this;
}

/* ----------------------------------------------------------------------
   returns whether this potential fn is radial or not
------------------------------------------------------------------------- */

int PotFns::is_radial() const
{
  return (fn_type_ == FnType::Radial);
}

/* ----------------------------------------------------------------------
   get maximum radial cutoff from all potentials
------------------------------------------------------------------------- */

double PotFns::get_max_rcut() const
{
  double max_rcut = 0.0;

  for (Basis *fn : fns)
    max_rcut = std::max(max_rcut, fn->get_max_rcut());

  return max_rcut;
}

/* ----------------------------------------------------------------------
   get maximum |y| value from all fns in this potential
------------------------------------------------------------------------- */

double PotFns::get_max_y_mag() const
{
  double max_y_mag = 0.0;

  for (Basis *fn : fns) {
    double fn_max_y_mag = fn->get_max_y_mag();
    if (std::abs( max_y_mag ) < std::abs( fn_max_y_mag ))
      max_y_mag = fn_max_y_mag;
  }

  return max_y_mag;
}

/* ----------------------------------------------------------------------
   get index of alloy basis fn knowing types of atoms
------------------------------------------------------------------------- */

int PotFns::get_alloy_idx(int typ_i, int typ_j) const
{
  switch(alloy_type_) {
  case AlloyType::Atom_i:
    return typ_i;
    break;
  case AlloyType::Atom_j:
    return typ_j;
    break;
  case AlloyType::Pair_ij:
    return (typ_i <= typ_j) ? typ_i * ntypes_ + typ_j - ((typ_i * (typ_i + 1)) / 2) :
                              typ_j * ntypes_ + typ_i - ((typ_j * (typ_j + 1)) / 2);
    break;
  default:
    return -1;
  }
}

/* ----------------------------------------------------------------------
   get the basis type common among all basis fns for this
   pot fn (return "undefined" if no common type)
------------------------------------------------------------------------- */

String PotFns::get_common_basis_type() const
{
  String common_basis_type = fns[0]->get_basis_type();

  for (Basis *fn : fns)
    if (common_basis_type != fn->get_basis_type())
      common_basis_type = "undefined";

  return common_basis_type;
}

/* ----------------------------------------------------------------------
   get number of adjustable coefficients in basis fns
------------------------------------------------------------------------- */

int PotFns::get_ncoeff() const
{
  int ncoeff = 0;

  for (Basis *fn : fns)
    ncoeff += fn->get_ncoeff();

  return ncoeff;
}

/* ----------------------------------------------------------------------
   return value of specified coeffient
------------------------------------------------------------------------- */

double& PotFns::coeff(int& idx)
{
  double *value = nullptr;
  for (Basis*& fn : fns) {
    value = &fn->coeff(idx);
    if (idx < 0) return *value;
  }

  return *value;
}

/* ----------------------------------------------------------------------
   perform modification to fn at specified coefficient
------------------------------------------------------------------------- */

void PotFns::modify_coeff(int idx, double height, double width)
{
  // Find basis function where idx lies in, and modify that basis fn
  int sum = 0;
  for (Basis*& fn : fns) {
    int fn_ncoeff = fn->get_ncoeff();
    sum += fn_ncoeff;
    if (idx < sum) {
      fn->modify_coeff(idx + fn_ncoeff - sum, height, width);
      return;
    }
  }

  return;
}

/* ----------------------------------------------------------------------
   check that a pot fn is similar enough to this
   pot fn for pairs/triplets to have identical info
------------------------------------------------------------------------- */

void PotFns::check_if_similar(PotFns& pot_fn) const
{
  // First check that each alloy fn is similar to corresponding alloy fn in other pot fn
  for (int i=0; i<n_fns; ++i)
    fns[i]->check_if_similar(pot_fn.fns[i]);  // error will be thrown if not similar

  // Then check that each alloy fn is similar to first alloy fn in other pot fn
  for (int j=1; j<n_fns; ++j)
    fns[0]->check_if_similar(pot_fn.fns[j]);  // error will be thrown if not similar

  return;
}

/* ----------------------------------------------------------------------
   refresh values in basis set to account for changes made earlier
------------------------------------------------------------------------- */

void PotFns::refresh_basis()
{
  for (Basis*& basis : fns)
    basis->refresh_basis();

  return;
}

/* ----------------------------------------------------------------------
   initialize size of alloy fns and fill with correct basis
   and dependency
      Atom_i,   // depends on atom in question, atom_i (A, B)
      Atom_j,   // depends on neighbor atom, atom_j (A, B)
      Pair_ij   // depends on pair of atoms i and j (AA, AB, BB)
------------------------------------------------------------------------- */

void PotFns::setup_pot_fns(int ntypes, AlloyType alloy_fn_type, FnType fn_type)
{
  alloy_type_ = alloy_fn_type;
  fn_type_ = fn_type;
  ntypes_ = ntypes;

  // Resize alloy_pots based on number of atom types and alloy fn type
  switch(alloy_type_) {
  case AlloyType::Atom_i:
    n_fns = ntypes;
    fns.resize(n_fns, nullptr);
    break;
  case AlloyType::Atom_j:
    n_fns = ntypes;
    fns.resize(n_fns, nullptr);
    break;
  case AlloyType::Pair_ij:
    n_fns = (ntypes * (ntypes + 1))/2;
    fns.resize(n_fns, nullptr);
    break;
  default:
    throw Error("ERROR: Ill-defined alloy type for potential function - ");
  }

  return;
}

/* ----------------------------------------------------------------------
   read in each alloy potential for this unique potential
------------------------------------------------------------------------- */

int PotFns::read_pot_fns(const StringList& lines, int line_num)
{
  for (Basis*& fn : fns) {
    if (lines[line_num].substr(0,2)!="#B") throw Error("ERROR: potential line (#B ...) missing - ");

    std::string label, basis;
    std::istringstream iss (lines[line_num++], std::istringstream::in);
    iss >> label >> basis;

    if (iss.fail())
      throw Error("ERROR: incorrect basis definition, use format (#B <basis type>) - ");

    // Set up basis type for this potential
    if (0) return line_num;
#define BASIS_CLASS
#define BasisStyle(key,Class) \
    else if (basis == #key) fn = new Class();
#include "style_basis.h"
#undef BasisStyle
#undef BASIS_CLASS
    else throw Error("ERROR: invalid basis type \"" + basis + "\" - ");

    line_num = fn->read_basis(lines, line_num);
  }

  return line_num;
}

/* ----------------------------------------------------------------------
   write out potential fn for each alloy fn
------------------------------------------------------------------------- */

void PotFns::write_pot_fns(std::ostream& oss) const
{
  for (Basis* fn : fns)
    fn->write_basis_fn(oss);

  return;
}

/* ----------------------------------------------------------------------
   write out potential fn for each alloy fn in lammps format
------------------------------------------------------------------------- */

void PotFns::write_lmp_fns(std::ostream& oss) const
{
  for (Basis* fn : fns)
    fn->write_lmp_basis_fn(oss);

  return;
}

/* ----------------------------------------------------------------------
   communicate data to all procs within group
------------------------------------------------------------------------- */

void PotFns::communicate(const Comm& comm, int root)
{
  for (Basis*& basis_fn : fns)
    basis_fn->communicate(comm, root);

  return;
}
