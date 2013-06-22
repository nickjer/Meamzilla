/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "spline3eq.h"
#include "error.h"
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Spline3Eq::Spline3Eq() : Spline(), step_(0.0), invstep_(0.0)
{
  basis_type_ = "spline3eq";

  // y' is only boundary condition for cubic splines
  grad_[0] = 0.0;  // left boundary condition
  is_grad_fixed_[0] = 1;
  grad_[1] = 0.0;  // right boundary condition
  is_grad_fixed_[1] = 1;
}

/* ---------------------------------------------------------------------- */

Spline3Eq::~Spline3Eq()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone spline
------------------------------------------------------------------------- */

Spline3Eq* Spline3Eq::clone() const
{
  return new Spline3Eq(*this);
}

/* ----------------------------------------------------------------------
   read in spline
------------------------------------------------------------------------- */

int Spline3Eq::read_basis(const StringList& lines, int line_num)
{
  // Read in spline definition for each function (<left> <right> <nknots>)
  std::istringstream iss (lines[line_num], std::istringstream::in);

  double left, right;
  int nknots;
  iss >> left >> right >> nknots;

  if (iss.fail())
    throw Error("ERROR: incorrect spline definition, use format (<left x> <right x> <nknots>) - ");

  set_spline(left, right, nknots);

  // Clear the string stream and reset it to next line
  iss.clear();
  iss.str(lines[++line_num]);

  // Get boundary condition constraints first: <left> <right>
  int l_constraint, r_constraint;
  iss >> l_constraint >> r_constraint;
  if (iss.fail())
    throw Error("ERROR: incorrect boundary condition constraint line, (<grad left fixed?> <grad right fixed?>) - ");
  if (l_constraint == -1) is_nat_[0] = 1;
  else if (l_constraint == 0) { is_grad_fixed_[0] = 0; ++ncoeff_; }
  else if (l_constraint == 1) is_grad_fixed_[0] = 1;
  else throw Error("ERROR: boundary condition constraint line has an incorrect constraint value for left side - ");
  if (r_constraint == -1) is_nat_[1] = 1;
  else if (r_constraint == 0) { is_grad_fixed_[1] = 0; ++ncoeff_; }
  else if (r_constraint == 1) is_grad_fixed_[1] = 1;
  else throw Error("ERROR: boundary condition constraint line has an incorrect constraint value for right side - ");

  // Get boundary conditions next: <left bc> <right bc>
  iss.clear();
  iss.str(lines[++line_num]);

  iss >> grad_[0] >> grad_[1];
  if (iss.fail())
    throw Error("ERROR: incorrect boundary condition line, format is (<left grad> <right grad>) - ");

  // Read in knots
  for (int k=0; k<nknots; ++k) {
    iss.clear();
    iss.str(lines[++line_num]);

    iss >> is_y_fixed_[k] >> y_[k];
    if (!is_y_fixed_[k]) ++ncoeff_;
    if (iss.fail())
      throw Error("ERROR: incorrect spline knot line, format is (<is y fixed?> <y-value>) - ");
  }

  return ++line_num;
}

/* ----------------------------------------------------------------------
   set new boundaries
------------------------------------------------------------------------- */

void Spline3Eq::readjust_x(double left, double right)
{
  // Reset 2nd derivs for interpolation later
  refresh_basis();

  Spline3Eq new_spline = *this;

  // Adjust step size
  new_spline.step_ = (right - left) / (nknots_ - 1);
  new_spline.invstep_ = 1.0 / new_spline.step_;

  // Set x-values
  new_spline.x_[0] = left;
  new_spline.x_[nknots_-1] = right;
  for (int i=1; i<nknots_-1; ++i)
    new_spline.x_[i] = left + i * new_spline.step_;

  // Set y-values
  for (int i=0; i<nknots_; ++i)
    new_spline.y_[i] = splint(new_spline.x_[i]);

  // Set y'-values
  if (!is_nat_[0])
    new_spline.grad_[0] = splint_grad(left);
  if (!is_nat_[1])
    new_spline.grad_[1] = splint_grad(right);

  new_spline.refresh_basis();

  *this = new_spline;

  return;
}

/* ----------------------------------------------------------------------
   rescale the x-coords by factor
------------------------------------------------------------------------- */

void Spline3Eq::rescale_x(double scale)
{
  step_ *= std::abs(scale);
  invstep_ = 1.0/step_;

  for (int i=0; i<nknots_; ++i)
    x_[i] *= scale;

  if (!is_nat_[0])
    grad_[0] /= scale;
  if (!is_nat_[1])
    grad_[1] /= scale;

  if (scale < 0.0) {
    std::reverse(x_.begin(),x_.end());
    std::reverse(y_.begin(),y_.end());
    std::reverse(&is_nat_[0],&is_nat_[1]);
    std::reverse(is_y_fixed_.begin(),is_y_fixed_.end());
    std::reverse(&grad_[0],&grad_[1]);
    std::reverse(&is_grad_fixed_[0],&is_grad_fixed_[1]);
  }
}

/* ----------------------------------------------------------------------
   Add linear fn a*x to spline
------------------------------------------------------------------------- */

void Spline3Eq::add_linear(double a)
{
  for (int i=0; i<nknots_; ++i)
    y_[i] += a * x_[i];

  if (!is_nat_[0])
    grad_[0] += a;
  if (!is_nat_[1])
    grad_[1] += a;
}

/* ----------------------------------------------------------------------
   refresh ypp values in spline to account for changes made earlier
------------------------------------------------------------------------- */

void Spline3Eq::refresh_basis()
{
  double p, qn, un;
  int n = nknots_;  // Number of spline knots
  std::vector<double> u(n);

  if (is_nat_[0]) {
    ypp_[0] = 0;
  } else {
    ypp_[0] = -0.5;
    u[0] = (3.0 * invstep_) * ((y_[1] - y_[0]) * invstep_ - grad_[0]);
  }

  for (int i = 1; i < n - 1; ++i) {
    p = 0.5 * ypp_[i - 1] + 2.0;
    ypp_[i] = (-0.5) / p;
    u[i] = (y_[i + 1] - y_[i]) * invstep_ - (y_[i] - y_[i - 1]) * invstep_;
    u[i] = (6.0 * u[i] * 0.5 * invstep_ - 0.5 * u[i - 1]) / p;
  }

  if (is_nat_[1]) {
    qn = un = 0.0;
  } else {
    qn = 0.5;
    un = (3.0 * invstep_) * (grad_[1] - (y_[n - 1] - y_[n - 2]) * invstep_);
  }

  ypp_[n - 1] = (un - qn * u[n - 2]) / (qn * ypp_[n - 2] + 1.0);

  for (int k = n - 2; k >= 0; --k)
    ypp_[k] = ypp_[k] * ypp_[k + 1] + u[k];

  return;
}

/* ----------------------------------------------------------------------
   communicate rest of spline data to other procs in group
------------------------------------------------------------------------- */

void Spline3Eq::communicate(const Comm& comm, int root)
{
  Spline::communicate(comm, root);

  comm.bcast(&step_, 1, MPI_DOUBLE, root);
  comm.bcast(&invstep_, 1, MPI_DOUBLE, root);
  comm.bcast(&ypp_[0], nknots_, MPI_DOUBLE, root);
  comm.bcast(&grad_[0], 2, MPI_DOUBLE, root);
  comm.bcast(&is_grad_fixed_[0], 2, MPI_INT, root);

  return;
}

/* ----------------------------------------------------------------------
   write out basis fn
------------------------------------------------------------------------- */

void Spline3Eq::write_basis_fn(std::ostream& oss) const
{
  oss << std::scientific;
  oss << std::setprecision(15);

  // Header for spline
  oss << "#B " << basis_type_ << std::endl;
  oss << x_[0] << " " << x_[nknots_-1] << " " << nknots_ << std::endl;

  // Output boundary condition constraints: <left> <right>
  // Followed by boundary conditions
  int constraint[2];
  for (int i=0; i<2; ++i)
    if (is_nat_[i]) constraint[i] = -1;
    else constraint[i] = is_grad_fixed_[i];
  oss << constraint[0] << " " << constraint[1] << std::endl;
  oss << grad_[0] << " " << grad_[1] << std::endl;

  // Output knot constraint and knot values
  for (int k=0; k<nknots_; ++k)
    oss << is_y_fixed_[k] << " " << y_[k] << std::endl;

  oss << std::endl;

  return;
}

/* ----------------------------------------------------------------------
   write out basis fn in lammps format
------------------------------------------------------------------------- */

void Spline3Eq::write_lmp_basis_fn(std::ostream& oss) const
{
  oss << std::scientific;
  oss << std::setprecision(15);

  // First line is number of knots
  oss << nknots_ << std::endl;

  // 2nd line is BC's
  oss << splint_grad(0, 0.0) << " " << splint_grad(nknots_-2, 1.0) << std::endl;

  // Output knot values
  for (int k=0; k<nknots_; ++k)
    oss << x_[k] << " " << y_[k] << " " << ypp_[k] << std::endl;

  return;
}

/* ----------------------------------------------------------------------
   compound assignment using multiplication with double
------------------------------------------------------------------------- */

Spline3Eq& Spline3Eq::operator*=(const double rhs)
{
  for (int i=0; i<nknots_; ++i)
    y_[i] *= rhs;

  if (!is_nat_[0])
    grad_[0] *= rhs;
  if (!is_nat_[1])
    grad_[1] *= rhs;

  return *this;
}

/* ----------------------------------------------------------------------
   compound assignment using division with double
------------------------------------------------------------------------- */

Spline3Eq& Spline3Eq::operator/=(const double rhs)
{
  for (int i=0; i<nknots_; ++i)
    y_[i] /= rhs;

  if (!is_nat_[0])
    grad_[0] /= rhs;
  if (!is_nat_[1])
    grad_[1] /= rhs;

  return *this;
}

/* ----------------------------------------------------------------------
   compound assignment using addition with another spline
------------------------------------------------------------------------- */

Spline3Eq& Spline3Eq::operator+=(const Basis& rhs)
{
  double rhs_x_max = rhs.get_max_rcut();

  for (int i=0; i<nknots_; ++i)
    y_[i] += (x_[i] < rhs_x_max) ? rhs.eval(x_[i]) : 0.;

  if (!is_nat_[0])
    grad_[0] += (x_[0] < rhs_x_max) ? rhs.eval_grad(x_[0]) : 0.;
  if (!is_nat_[1])
    grad_[1] += (x_.back() < rhs_x_max) ? rhs.eval_grad(x_.back()) : 0.;

  return *this;
}

/* ----------------------------------------------------------------------
   setup equidistant spline
------------------------------------------------------------------------- */

void Spline3Eq::set_spline(double left, double right, int nknots)
{
  nknots_ = nknots;
  x_.resize(nknots);
  y_.resize(nknots);
  ypp_.resize(nknots);
  is_y_fixed_.resize(nknots,0);

  step_ = (right - left)/(nknots - 1);
  invstep_ = 1./step_;

  x_[0] = left;
  x_[nknots-1] = right;

  for (int i=1; i<nknots-1; ++i)
    x_[i] = left + step_*i;
}

/* ----------------------------------------------------------------------
   return BC if it is specified coefficient
------------------------------------------------------------------------- */

double& Spline3Eq::bc_coeff(int& idx)
{
  for (int i=0; i<2; ++i)
    if (!is_grad_fixed_[i] && !is_nat_[i]) if (idx-- == 0) return grad_[i];

  return y_[0];
}

/* ----------------------------------------------------------------------
   perform modification to basis at specified coefficient
------------------------------------------------------------------------- */

void Spline3Eq::modify_bc_coeff(int idx, double height)
{
  // Check gradients and just add the height, then leave
  for (int i=0; i<2; ++i)
    if (!is_grad_fixed_[i] && !is_nat_[i]) if (idx-- == 0) {
      grad_[i] += height;
      return;
    }

  return;
}
