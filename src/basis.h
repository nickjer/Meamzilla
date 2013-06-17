#ifndef MEAMZ_BASIS_H
#define MEAMZ_BASIS_H

#include "comm.h"

#include <string>
#include <vector>

namespace MEAMZ_NS
{

typedef std::string String;
typedef std::vector<String> StringList;

class Basis
{
public:
  Basis();
  virtual ~Basis();

  virtual Basis* clone() const = 0;                     // Copy polymorphic objects

  virtual double get_min_rcut() const = 0;              // Get minimum radius for this basis fn
  virtual double get_max_rcut() const = 0;              // Get maximum radius for this basis fn
  String get_basis_type() const;                        // Get basis type for this basis fn
  virtual double get_max_y_mag() const = 0;             // Get maximum |y| value for this basis fn
  virtual int read_basis(const StringList&, int) = 0;   // Read in basis fn from StringList
  int get_ncoeff() const;                               // Get number of adjustable coefficients for this basis fn
  virtual double& coeff(int&) = 0;                      // Return value of specified coeffient
  virtual void modify_coeff(int, double, double) = 0;   // Perform modification to basis at specified coefficient

  virtual double eval(double) const = 0;                // Evaluate the basis fn at defined point
  virtual double eval_grad(double) const = 0;           // Evaluate gradient of basis fn at defined point
  virtual double eval_comb(double, double *) const = 0; // Evaluate basis fn and grad together
  virtual double eval_grad2(double) const = 0;          // Evaluate 2nd deriv of basis fn at defined point

  virtual void readjust_x(double, double) = 0;          // Set new boundaries
  virtual void rescale_x(double) = 0;                   // Rescale the x-coords by factor
  virtual void add_linear(double) = 0;                  // Add linear fn a*x to spline

  virtual void check_if_similar(Basis *) const = 0;     // Check that a basis fn is similar enough to this
                                                        // basis fn for pairs/triplets to have identical info

  virtual void refresh_basis() = 0;                     // Refresh values in basis set to account for changes made earlier
  virtual void communicate(const Comm&, int) = 0;       // Communicate basis fn to other procs in group

  virtual void write_basis_fn(std::ostream&) const = 0;     // Write out basis fn
  virtual void write_lmp_basis_fn(std::ostream&) const = 0; // Write out basis fn in lammps format

  virtual Basis& operator*=(const double) = 0;          // Compound assignment using multiplication with double
  virtual Basis& operator/=(const double) = 0;          // Compound assignment using division with double
  virtual Basis& operator+=(const Basis&) = 0;          // Compound assignment using addition with another spline

protected:
  int ncoeff_;          // Number of adjustable coefficients that are not fixed

  String basis_type_;   // Basis type name

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_BASIS_H
