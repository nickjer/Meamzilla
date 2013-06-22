/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef BASIS_CLASS

BasisStyle(spline3eq,Spline3Eq)

#else

#ifndef MEAMZ_SPLINE3EQ_H
#define MEAMZ_SPLINE3EQ_H

#include "spline.h"

namespace MEAMZ_NS
{

class Spline3Eq : public Spline
{
public:
  Spline3Eq();
  virtual ~Spline3Eq();

  virtual Spline3Eq* clone() const;                       // Copy polymorphic objects

  virtual int read_basis(const StringList&, int);         // Read in spline from StringList

  virtual void readjust_x(double, double);                // Set new boundaries
  virtual void rescale_x(double);                         // Rescale the x-coords by factor
  virtual void add_linear(double);                        // Add linear fn a*x to spline

  virtual void refresh_basis();                           // Refresh ypp values in spline to account for changes made earlier
  virtual void communicate(const Comm&, int);             // Communicate basis fn to other procs in group

  virtual void write_basis_fn(std::ostream&) const;       // Write out basis fn
  virtual void write_lmp_basis_fn(std::ostream&) const;   // Write out basis fn

  double splint(int, double) const;                       // Interpolate fn from splines given knot_idx & shift
  double splint_grad(int, double) const;                  // Interpolate grad from splines given knot_idx & shift
  double splint_comb(int, double, double *) const;        // Interpolate fn & grad from splines given knot_idx & shift
  double splint_grad2(int, double) const;                 // Interpolate 2nd deriv from splines given knot_idx & shift

  double splint(double) const;                            // Interpolate fn from point
  double splint_grad(double) const;                       // Interpolate grad from point
  double splint_comb(double, double *) const;             // Interpolate fn & grad from point
  double splint_grad2(double) const;                      // Interpolate 2nd deriv from point

  virtual Spline3Eq& operator*=(const double);            // Compound assignment using multiplication with double
  virtual Spline3Eq& operator/=(const double);            // Compound assignment using division with double
  virtual Spline3Eq& operator+=(const Basis&);            // Compound assignment using addition with another spline

protected:
  double step_;                 // distance between spline knots
  double invstep_;              // inverse of spline knot distance

  std::vector<double> ypp_;     // y'' at each knot

  double grad_[2];              // gradient at the end points
  int is_grad_fixed_[2];        // Whether gradient is fixed or not

  virtual void set_spline(double, double, int);     // Set up equidistant spline
  virtual double& bc_coeff(int&);                   // Return BC if it is specified coefficient
  virtual void modify_bc_coeff(int, double);        // Perform modification to boundary condition

private:

};

/* ----------------------------------------------------------------------
   Interpolate fn from splines given knot_idx & shift
------------------------------------------------------------------------- */
inline
double Spline3Eq::splint(int knot_idx, double shift) const
{
  double a = 1.0 - shift;
  double p1 = y_[knot_idx];
  double d21 = ypp_[knot_idx];
  double p2 = y_[++knot_idx];
  double d22 = ypp_[knot_idx];

  return a*p1 + shift*p2 + ((a*a*a-a)*d21 + (shift*shift*shift-shift)*d22)*step_*step_/6.0;
}

/* ----------------------------------------------------------------------
   Interpolate grad from splines given knot_idx & shift
------------------------------------------------------------------------- */
inline
double Spline3Eq::splint_grad(int knot_idx, double shift) const
{
  double a = 1.0 - shift;
  double p1 = y_[knot_idx];
  double d21 = ypp_[knot_idx];
  double p2 = y_[++knot_idx];
  double d22 = ypp_[knot_idx];

  return (p2-p1)*invstep_ + ((3*shift*shift-1)*d22 - (3*a*a-1)*d21)*step_/6.0;
}

/* ----------------------------------------------------------------------
   Interpolate fn & grad from splines given knot_idx & shift
------------------------------------------------------------------------- */
inline
double Spline3Eq::splint_comb(int knot_idx, double shift, double *grad) const
{
  double a = 1.0 - shift;
  double p1 = y_[knot_idx];
  double d21 = ypp_[knot_idx];
  double p2 = y_[++knot_idx];
  double d22 = ypp_[knot_idx];

  *grad = (p2-p1)*invstep_ + ((3*shift*shift-1)*d22 - (3*a*a-1)*d21)*step_/6.0;

  return a*p1 + shift*p2 + ((a*a*a-a)*d21 + (shift*shift*shift-shift)*d22)*step_*step_/6.0;
}

/* ----------------------------------------------------------------------
   Interpolate 2nd deriv from splines given knot_idx & shift
------------------------------------------------------------------------- */
inline
double Spline3Eq::splint_grad2(int knot_idx, double shift) const
{
  double a = 1.0 - shift;
  double d21 = ypp_[knot_idx];
  double d22 = ypp_[++knot_idx];

  return a*d21 + shift*d22;
}

/* ----------------------------------------------------------------------
   Interpolate fn from point
------------------------------------------------------------------------- */
inline
double Spline3Eq::splint(double r) const
{
  if (r < x_[0])
    return y_[0] + grad_[0] * (r - x_[0]);
  else if (r > x_[nknots_-1])
    return y_[nknots_-1] + grad_[1] * (r - x_[nknots_-1]);

  int k;
  double delta_r = r - x_[0];

  // Stay within spline lower bound
  if (delta_r < 0) {
    k = 0;
  } else {
    k = int(delta_r * invstep_);
    if ( k >= nknots_-1 ) k = nknots_-2;  // stay within one less than spline upper bound
  }

  double b = (delta_r - k*step_)*invstep_;
  double a = 1.0 - b;
  double p1 = y_[k];
  double d21 = ypp_[k];
  double p2 = y_[++k];
  double d22 = ypp_[k];

  return a*p1 + b*p2 + ((a*a*a-a)*d21 + (b*b*b-b)*d22)*(step_*step_)/6.0;
}

/* ----------------------------------------------------------------------
   Interpolate grad from point
------------------------------------------------------------------------- */
inline
double Spline3Eq::splint_grad(double r) const
{
  if (r < x_[0])
    return grad_[0];
  else if (r > x_[nknots_-1])
    return grad_[1];

  int k;
  double delta_r = r - x_[0];

  // Stay within spline lower bound
  if (delta_r < 0) {
    k = 0;
  } else {
    k = int(delta_r * invstep_);
    if ( k >= nknots_-1 ) k = nknots_-2;  // stay within one less than spline upper bound
  }

  double b = (delta_r - k*step_)*invstep_;
  double a = 1.0 - b;
  double p1 = y_[k];
  double d21 = ypp_[k];
  double p2 = y_[++k];
  double d22 = ypp_[k];

  return (p2-p1)*invstep_ + ((3*(b*b)-1)*d22 - (3*(a*a)-1)*d21) * step_/6.0;
}

/* ----------------------------------------------------------------------
   Interpolate fn & grad from point
------------------------------------------------------------------------- */
inline
double Spline3Eq::splint_comb(double r, double *grad) const
{
  if (r < x_[0]) {
    *grad = grad_[0];
    return y_[0] + grad_[0] * (r - x_[0]);
  } else if (r > x_[nknots_-1]) {
    *grad = grad_[1];
    return y_[nknots_-1] + grad_[1] * (r - x_[nknots_-1]);
  }

  int k;
  double delta_r = r - x_[0];

  // Stay within spline lower bound
  if (delta_r < 0) {
    k = 0;
  } else {
    k = int(delta_r * invstep_);
    if ( k >= nknots_-1 ) k = nknots_-2;  // stay within one less than spline upper bound
  }

  double b = (delta_r - k*step_)*invstep_;
  double a = 1.0 - b;
  double p1 = y_[k];
  double d21 = ypp_[k];
  double p2 = y_[++k];
  double d22 = ypp_[k];

  *grad = (p2-p1)*invstep_ + ((3*(b*b)-1)*d22 - (3*(a*a)-1)*d21) * step_/6.0;

  return a*p1 + b*p2 + ((a*a*a-a)*d21 + (b*b*b-b)*d22)*(step_*step_)/6.0;
}

/* ----------------------------------------------------------------------
   Interpolate 2nd deriv from point
------------------------------------------------------------------------- */
inline
double Spline3Eq::splint_grad2(double r) const
{
  if (r < x_[0])
    return 0.0;
  else if (r > x_[nknots_-1])
    return 0.0;

  int k;
  double delta_r = r - x_[0];

  // Stay within spline lower bound
  if (delta_r < 0) {
    k = 0;
  } else {
    k = int(delta_r * invstep_);
    if ( k >= nknots_-1 ) k = nknots_-2;  // stay within one less than spline upper bound
  }

  double b = (delta_r - k*step_)*invstep_;
  double a = 1.0 - b;
  double d21 = ypp_[k];
  double d22 = ypp_[++k];

  return a*d21 + b*d22;
}

} // namespace MEAMZ_NS

#endif // MEAMZ_SPLINE3EQ_H
#endif // BASIS_CLASS
