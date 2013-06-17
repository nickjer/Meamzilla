#ifndef MEAMZ_SPLINE_H
#define MEAMZ_SPLINE_H

#include "basis.h"
#include <vector>

namespace MEAMZ_NS
{

class Spline : public Basis
{
public:
  Spline();
  virtual ~Spline();

  virtual Spline* clone() const = 0;                        // Copy polymorphic objects

  double get_min_rcut() const;                              // Get minimum radius for this spline (first x-knot)
  double get_max_rcut() const;                              // Get maximum radius for this spline (last x-knot)
  double get_max_y_mag() const;                             // Get maximum |y| value for this spline
  virtual int read_basis(const StringList&, int) = 0;       // Read in spline from StringList
  int get_knotshift(double, double&) const;                 // Get knot and shift of x-coord 'r'

  double& coeff(int&);                                      // Return value of specified coeffient
  void modify_coeff(int, double, double);                   // Perform modification to basis at specified coefficient

  double eval(double) const;                                // Evaluate the basis fn at defined point
  double eval_grad(double) const;                           // Evaluate gradient of basis fn at defined point
  double eval_comb(double, double *) const;                 // Evaluate basis fn and grad together
  double eval_grad2(double) const;                          // Evaluate 2nd deriv of basis fn at defined point

  virtual void readjust_x(double, double) = 0;              // Set new boundaries
  virtual void rescale_x(double) = 0;                       // Rescale the x-coords by factor
  virtual void add_linear(double) = 0;                      // Add linear fn a*x to spline

  void check_if_similar(Basis *) const;                     // Check that a basis fn is similar enough to this
                                                            // basis fn for pairs/triplets to have identical info
  virtual void refresh_basis() = 0;                         // Refresh values in spline to account for changes made earlier
  virtual void communicate(const Comm&, int);               // Communicate basis fn to other procs in group

  virtual void write_basis_fn(std::ostream&) const = 0;     // Write out basis fn
  virtual void write_lmp_basis_fn(std::ostream&) const = 0; // Write out basis fn in lammps format

  virtual double splint(int, double) const;                 // Interpolate fn from splines given knot_idx & shift
  virtual double splint_grad(int, double) const;            // Interpolate grad from splines given knot_idx & shift
  virtual double splint_comb(int, double, double *) const;  // Interpolate fn & grad from splines given knot_idx & shift

  virtual double splint(double) const;                      // Interpolate fn from point
  virtual double splint_grad(double) const;                 // Interpolate grad from point
  virtual double splint_comb(double, double *) const;       // Interpolate fn & grad from point
  virtual double splint_grad2(double) const;                // Interpolate 2nd deriv from point

  virtual Spline& operator*=(const double) = 0;             // Compound assignment using multiplication with double
  virtual Spline& operator/=(const double) = 0;             // Compound assignment using division with double
  virtual Spline& operator+=(const Basis&) = 0;             // Compound assignment using addition with another spline

protected:
  int nknots_;                        // size of potential
  std::vector<double> x_, y_;         // x and y values at each spline knot

  int is_nat_[2];                     // is natural boundary?
  std::vector<int> is_y_fixed_;       // Whether y-value is fixed (at each knot)

  virtual double& bc_coeff(int&) = 0;             // Return BC if it is specified coefficient
  virtual void modify_bc_coeff(int, double) = 0;  // Perform modification to boundary condition

private:

};

/* ----------------------------------------------------------------------
   Interpolate fn from splines given knot_idx & shift
------------------------------------------------------------------------- */
inline
double Spline::splint(int knot_idx, double shift) const
{
  return splint(knot_idx, shift);
}

/* ----------------------------------------------------------------------
   Interpolate grad from splines given knot_idx & shift
------------------------------------------------------------------------- */
inline
double Spline::splint_grad(int knot_idx, double shift) const
{
  return splint_grad(knot_idx, shift);
}

/* ----------------------------------------------------------------------
   Interpolate fn & grad from splines given knot_idx & shift
------------------------------------------------------------------------- */
inline
double Spline::splint_comb(int knot_idx, double shift, double *grad) const
{
  return splint_comb(knot_idx, shift, grad);
}

/* ----------------------------------------------------------------------
   Interpolate fn from point
------------------------------------------------------------------------- */
inline
double Spline::splint(double r) const
{
  return splint(r);
}

/* ----------------------------------------------------------------------
   Interpolate grad from point
------------------------------------------------------------------------- */
inline
double Spline::splint_grad(double r) const
{
  return splint_grad(r);
}

/* ----------------------------------------------------------------------
   Interpolate fn & grad from point
------------------------------------------------------------------------- */
inline
double Spline::splint_comb(double r, double *grad) const
{
  return splint_comb(r, grad);
}

/* ----------------------------------------------------------------------
   Interpolate 2nd deriv from point
------------------------------------------------------------------------- */
inline
double Spline::splint_grad2(double r) const
{
  return splint_grad2(r);
}

} // namespace MEAMZ_NS

#endif // MEAMZ_SPLINE_H
