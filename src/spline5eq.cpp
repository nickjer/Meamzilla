/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "spline5eq.h"
#include "matrix.h"
#include "error.h"
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

extern "C" {
void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
}

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Spline5Eq::Spline5Eq() : Spline3Eq()
{
  basis_type_ = "spline5eq";

  // y'' is a boundary condition for quintic splines
  grad2_[0] = 0.0;  // left boundary condition
  grad2_[1] = 0.0;  // right boundary condition
}

/* ---------------------------------------------------------------------- */

Spline5Eq::~Spline5Eq()
{
  //dtor
}

/* ----------------------------------------------------------------------
   clone spline
------------------------------------------------------------------------- */

Spline5Eq* Spline5Eq::clone() const
{
  return new Spline5Eq(*this);
}

/* ----------------------------------------------------------------------
   read in spline
------------------------------------------------------------------------- */

int Spline5Eq::read_basis(const StringList& lines, int line_num)
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

  iss.clear();
  iss.str(lines[++line_num]);

  iss >> grad2_[0] >> grad2_[1];
  if (iss.fail())
    throw Error("ERROR: incorrect boundary condition line, format is (<left grad2> <right grad2>) - ");

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

void Spline5Eq::readjust_x(double left, double right)
{
  // Reset 2nd derivs for interpolation later
  refresh_basis();

  Spline5Eq new_spline = *this;

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

  // Set y' & y''-values
  if (!is_nat_[0]) {
    new_spline.grad_[0] = splint_grad(left);
    new_spline.grad2_[0] = splint_grad2(left);
  }
  if (!is_nat_[1]) {
    new_spline.grad_[1] = splint_grad(right);
    new_spline.grad2_[1] = splint_grad2(right);
  }

  new_spline.refresh_basis();

  *this = new_spline;

  return;
}

/* ----------------------------------------------------------------------
   rescale the x-coords by factor
------------------------------------------------------------------------- */

void Spline5Eq::rescale_x(double scale)
{
  step_ *= std::abs(scale);
  invstep_ = 1.0/step_;

  for (int i=0; i<nknots_; ++i)
    x_[i] *= scale;

  if (!is_nat_[0]) {
    grad_[0] /= scale;
    grad2_[0] /= scale*scale;
  }
  if (!is_nat_[1]) {
    grad_[1] /= scale;
    grad2_[1] /= scale*scale;
  }

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

void Spline5Eq::add_linear(double a)
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

void Spline5Eq::refresh_basis()
{
  int n = nknots_;
  int ndim = 2*n;
  Matrix<double> A(ndim,ndim,0.0); // A[row][col]
  std::vector<double> b(ndim);     // b[row][1] <- vector

  // Set the meat of the 4 block matrices
  for (int i=1; i<n-1; ++i) { // rows
    A(i,i-1) =  1./12;
    A(i,i)   = 10./12;
    A(i,i+1) =  1./12;

    A(i,n+i-1) =  1./180;
    A(i,n+i)   = -1./90;
    A(i,n+i+1) =  1./180;

    A(n+i,i-1) =  1.;
    A(n+i,i)   = -2.;
    A(n+i,i+1) =  1.;

    A(n+i,n+i-1) = 1./6;
    A(n+i,n+i)   = 2./3;
    A(n+i,n+i+1) = 1./6;

    // Meat of b-vector
    b[i]   = y_[i+1]-2*y_[i]+y_[i-1];
    b[n+i] = 0.;
  }

  // Next set BC's
  if (is_nat_[0]) {
    // S0''' = 0
    A(0,0) = -1.;
    A(0,1) =  1.;
    A(0,n)   = 1./3;
    A(0,n+1) = 1./6;
    b[0] = 0.;

    // S0'''' = 0
    A(n,n) = 1.;
    b[n] = 0.;
  } else {
    // S0' = grad[0]
    A(0,0) = 1./3;
    A(0,1) = 1./6;
    A(0,n)   = 1./45;
    A(0,n+1) = 7./360;
    b[0] = y_[1]-y_[0]-grad_[0]*step_;

    // S0'' = grad2[0]
    A(n,0) = 1.;
    b[n] = grad2_[0]*(step_*step_);
  }

  if (is_nat_[1]) {
    // SN'''' = 0
    A(n-1,n-2) = -1.;
    A(n-1,n-1) =  1.;
    A(n-1,2*n-2) = -1./6;
    A(n-1,2*n-1) = -1./3;
    b[n-1] = 0.;

    // SN'''' = 0
    A(2*n-1,2*n-1) = 1.;
    b[2*n-1] = 0.;
  } else {
    // SN' = grad[1]
    A(n-1,n-2) = -1./6;
    A(n-1,n-1) = -1./3;
    A(n-1,2*n-2) = -7./360;
    A(n-1,2*n-1) = -1./45;
    b[n-1] = y_[n-1]-y_[n-2]-grad_[1]*step_;

    // SN'' = grad2[1]
    A(2*n-1,n-1) = 1.;
    b[2*n-1] = grad2_[1]*(step_*step_);
  }

  // Solve for A.x = b
  std::vector<int> ipiv(ndim); // Keeps track of LU pivoting
  int info;                    // Info integer taken from matrix solver
  int nrhs = 1;                // The number of right hand sides, i.e., the number of columns of the matrix B
  dgesv_(&ndim, &nrhs, &A(0,0), &ndim, &ipiv[0], &b[0], &ndim, &info);

/*
  if (info > 0 && info <= n) {
    std::ostringstream oss;
    oss << "Linear equation system A is singular because U[i][i] = 0 at i = " << info;
    if (mmz->universe->me == 0) mmz->error->all(oss.str());
  }
*/

  for (int i=0; i<n; ++i) {
    ypp_[i] = b[i]*(invstep_*invstep_);
    yp4_[i] = -b[n+i]*(invstep_*invstep_*invstep_*invstep_);
  }

  return;
}

/* ----------------------------------------------------------------------
   communicate rest of spline data to other procs in group
------------------------------------------------------------------------- */

void Spline5Eq::communicate(const Comm& comm, int root)
{
  Spline3Eq::communicate(comm, root);

  comm.bcast(&yp4_[0], nknots_, MPI_DOUBLE, root);
  comm.bcast(&grad2_[0], 2, MPI_DOUBLE, root);

  return;
}

/* ----------------------------------------------------------------------
   write out basis fn
------------------------------------------------------------------------- */

void Spline5Eq::write_basis_fn(std::ostream& oss) const
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
  oss << grad2_[0] << " " << grad2_[1] << std::endl;

  // Output knot constraint and knot values
  for (int k=0; k<nknots_; ++k)
    oss << is_y_fixed_[k] << " " << y_[k] << std::endl;

  oss << std::endl;

  return;
}

/* ----------------------------------------------------------------------
   write out basis fn in lammps format
------------------------------------------------------------------------- */

void Spline5Eq::write_lmp_basis_fn(std::ostream& oss) const
{
  oss << std::scientific;
  oss << std::setprecision(15);

  // First line is number of knots
  oss << nknots_ << std::endl;

  // 2nd line is BC's
  oss << splint_grad(0, 0.0) << " " << splint_grad(nknots_-2, 1.0) << std::endl;

  // 3rd line is meaningless
  oss << "0 0 0 0" << std::endl;

  // Output knot values
  for (int k=0; k<nknots_; ++k)
    oss << x_[k] << " " << y_[k] << " " << ypp_[k] << " " << yp4_[k] << std::endl;

  return;
}

/* ----------------------------------------------------------------------
   compound assignment using multiplication with double
------------------------------------------------------------------------- */

Spline5Eq& Spline5Eq::operator*=(const double rhs)
{
  for (int i=0; i<nknots_; ++i)
    y_[i] *= rhs;

  if (!is_nat_[0]) {
    grad_[0] *= rhs;
    grad2_[0] *= rhs;
  }
  if (!is_nat_[1]) {
    grad_[1] *= rhs;
    grad2_[1] *= rhs;
  }

  return *this;
}

/* ----------------------------------------------------------------------
   compound assignment using division with double
------------------------------------------------------------------------- */

Spline5Eq& Spline5Eq::operator/=(const double rhs)
{
  for (int i=0; i<nknots_; ++i)
    y_[i] /= rhs;

  if (!is_nat_[0]) {
    grad_[0] /= rhs;
    grad2_[0] /= rhs;
  }
  if (!is_nat_[1]) {
    grad_[1] /= rhs;
    grad2_[1] /= rhs;
  }

  return *this;
}

/* ----------------------------------------------------------------------
   compound assignment using addition with another spline
------------------------------------------------------------------------- */

Spline5Eq& Spline5Eq::operator+=(const Basis& rhs)
{
  double rhs_x_max = rhs.get_max_rcut();

  for (int i=0; i<nknots_; ++i)
    y_[i] += (x_[i] < rhs_x_max) ? rhs.eval(x_[i]) : 0.;

  if (!is_nat_[0]) {
    grad_[0]  += (x_[0] < rhs_x_max) ? rhs.eval_grad(x_[0]) : 0.;
    grad2_[0] += (x_[0] < rhs_x_max) ? rhs.eval_grad2(x_[0]) : 0.;
  }
  if (!is_nat_[1]) {
    grad_[1]  += (x_.back() < rhs_x_max) ? rhs.eval_grad(x_.back()) : 0.;
    grad2_[1] += (x_.back() < rhs_x_max) ? rhs.eval_grad2(x_.back()) : 0.;
  }

  return *this;
}

/* ----------------------------------------------------------------------
   setup equidistant spline
------------------------------------------------------------------------- */

void Spline5Eq::set_spline(double left, double right, int nknots)
{
  Spline3Eq::set_spline(left, right, nknots);

  yp4_.resize(nknots);

  return;
}

/* ----------------------------------------------------------------------
   return BC if it is specified coefficient
------------------------------------------------------------------------- */

double& Spline5Eq::bc_coeff(int& idx)
{
  for (int i=0; i<2; ++i) {
    if (!is_grad_fixed_[i] && !is_nat_[i]) {
      if (idx-- == 0) return grad_[i];
      if (idx-- == 0) return grad2_[i];
    }
  }

  return y_[0];
}

/* ----------------------------------------------------------------------
   perform modification to basis at specified coefficient
------------------------------------------------------------------------- */

void Spline5Eq::modify_bc_coeff(int idx, double height)
{
  // Check gradients and just add the height, then leave
  for (int i=0; i<2; ++i)
    if (!is_grad_fixed_[i] && !is_nat_[i]) {
      if (idx-- == 0) {
        grad_[i] += height;
        return;
      }
      if (idx-- == 0) {
        grad2_[i] += height;
        return;
      }
    }

  return;
}
