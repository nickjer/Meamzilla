/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "spline.h"
#include "error.h"

#include <cmath>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Spline
::Spline() : Basis(), nknots_(0)
{
  // Initialize whether natural BC's on left and right
  is_nat_[0] = 0;
  is_nat_[1] = 0;
}

/* ---------------------------------------------------------------------- */

Spline::~Spline()
{
  //dtor
}

/* ----------------------------------------------------------------------
   get minimum radius for this spline (first x-knot)
------------------------------------------------------------------------- */

double Spline::get_min_rcut() const
{
  return x_[0];
}

/* ----------------------------------------------------------------------
   get maximum radius for this spline (last x-knot)
------------------------------------------------------------------------- */

double Spline::get_max_rcut() const
{
  return x_.back();
}

/* ----------------------------------------------------------------------
   get maximum |y| value for this spline
------------------------------------------------------------------------- */

double Spline::get_max_y_mag() const
{
  double max_y_mag = 0.0;

  for (double y : y_)
    if (std::abs(max_y_mag) < std::abs(y))
      max_y_mag = y;

  return max_y_mag;
}

/* ----------------------------------------------------------------------
   get lower bound knot and
   shift (percentage shift from lower bound knot) of x-coord 'r'
------------------------------------------------------------------------- */

int Spline::get_knotshift(double r, double& shift) const
{
  // Check if pair is inside pair potential first
  // then compute its position along spline
  int knot = -1;
  if (r < x_.back() * (1. + 1.e-8)) { // due to slight numerical inaccuracies (especially with cos(theta))
    int klo = 0;
    int khi = nknots_ - 1;

    // Find index by bisection
    while (khi - klo > 1) {
      int k = (khi + klo) >> 1;
      if (x_[k] > r)
        khi = k;
      else
        klo = k;
    }

    knot = klo; // set the lower bound knot value

    double h = x_[khi] - x_[klo];
    shift = (r - x_[klo]) / h;  // set the percentage shift from lower bound knot
  } // if pair is inside radial cutoff

  return knot;
}

/* ----------------------------------------------------------------------
   return value of specified coeffient
------------------------------------------------------------------------- */
#include <iostream>
double& Spline::coeff(int& idx)
{
  // First check through each spline knot
  for (int k=0; k<nknots_; ++k)
    if (!is_y_fixed_[k]) if (idx-- == 0) return y_[k];

  // Next check boundary conditions
  double& value = bc_coeff(idx);
  return value;
}

/* ----------------------------------------------------------------------
   perform modification to basis at specified coefficient
------------------------------------------------------------------------- */

void Spline::modify_coeff(int idx, double height, double width)
{
  // Find coefficient first
  int found = 0;
  int idx_k = -1;
  for (int k=0; k<nknots_; ++k) {
    if (!is_y_fixed_[k]) if (idx-- == 0) {
      idx_k = k;
      found = 1;
      break;
    }
  }

  // User wants to modify a BC instead of a knot
  if (found == 0) {
    modify_bc_coeff(idx, height);
    return;
  }

  // Found a spline knot to modify
  const double INV_SQRT_2PI = 0.3989422804014326779;  // used in Gaussian distribution

  for (int i=-int(4*width); i<=int(4*width); ++i) {
    int k = idx_k + i;
    if (k >= 0 && k < nknots_)
      if (!is_y_fixed_[k])
        y_[k] += INV_SQRT_2PI*std::exp(-i*i*0.5/(width*width)) * height;  // Gaussian distribution bump
  }
  return;
}

/* ----------------------------------------------------------------------
   evaluate the basis fn at defined point
------------------------------------------------------------------------- */

double Spline::eval(double r) const
{
  return splint(r);
}

/* ----------------------------------------------------------------------
   evaluate gradient of basis fn at defined point
------------------------------------------------------------------------- */

double Spline::eval_grad(double r) const
{
  return splint_grad(r);
}

/* ----------------------------------------------------------------------
   evaluate basis fn and grad together
------------------------------------------------------------------------- */

double Spline::eval_comb(double r, double *grad) const
{
  return splint_comb(r, grad);
}

/* ----------------------------------------------------------------------
   evaluate 2nd deriv of basis fn at defined point
------------------------------------------------------------------------- */

double Spline::eval_grad2(double r) const
{
  return splint_grad2(r);
}

/* ----------------------------------------------------------------------
   check that a basis fn is similar enough to this
   basis fn for pairs/triplets to have identical info
------------------------------------------------------------------------- */

void Spline::check_if_similar(Basis *basis_ptr) const
{
  Spline& basis = *static_cast<Spline*>(basis_ptr);

  if (nknots_ != basis.nknots_)
    throw Error("ERROR: Number of knots don't match - ");

  for (int i=0; i<nknots_; ++i)
    if (x_[i] != basis.x_[i])
      throw Error("ERROR: Spacing of knots don't match - ");

  return;
}

/* ----------------------------------------------------------------------
   communicate basis fn to other procs in group
------------------------------------------------------------------------- */

void Spline::communicate(const Comm& comm, int root)
{
  comm.bcast(&nknots_, 1, MPI_INT, root);
  if (comm.get_rank() != root) {
    x_.resize(nknots_);
    y_.resize(nknots_);
    is_y_fixed_.resize(nknots_);
  }
  comm.bcast(&x_[0], nknots_, MPI_DOUBLE, root);
  comm.bcast(&y_[0], nknots_, MPI_DOUBLE, root);
  comm.bcast(&is_y_fixed_[0], nknots_, MPI_INT, root);
  comm.bcast(&is_nat_[0], 2, MPI_INT, root);

  return;
}
