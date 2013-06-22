/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef OPTIMIZER_CLASS

OptimizerStyle(powell,OptPowell)

#else

#ifndef MEAMZ_OPT_POWELL_H
#define MEAMZ_OPT_POWELL_H

#include "optimizer.h"
#include "matrix.h"
#include "line_min_brent.h"

namespace MEAMZ_NS
{

class OptPowell : public Optimizer
{
public:
  OptPowell(class Meamzilla *, std::ostream * = nullptr);
  virtual ~OptPowell();

  virtual void init();              // Initialize the optimizer
  virtual void compute(PotList&);   // Use the optimizer on a potential list

  virtual Potential* compute_pot(const Potential *, int) const;  // Optimization loop

protected:
  double d_eps_;  // Error margin for powell minimization; if improvement is smaller, optimization is stopped

  // Constants used in the Powell minimization algorithm
  static constexpr double EPS = 1.E-4;
  static constexpr double PRECISION = 1.E-7;
  static constexpr double NOTHING = 1.E-12;
  static constexpr int INNERLOOPS = 801;
  static constexpr double TOOBIG = 1.E+4;

  void print_powell(int, double, int, int, int) const;      // Print results

  int gamma_init(const Potential *, const ErrorVec&, int, int, Matrix<double> *, Matrix<double> *) const;
  int gamma_update(int, Potential *, const ErrorVec&, int, int, LineMinBrent&, Vector<double> *, Matrix<double> *) const;
  void lineqsys_init(const ErrorVec&, int, int, const Matrix<double>&, Matrix<double> *, Vector<double> *) const;
  void lineqsys_update(int, const ErrorVec&, double, double, const Matrix<double>&, Matrix<double> *, Vector<double> *) const;

private:

};

} // MEAMZ_NS

#endif // MEAMZ_OPT_POWELL_H
#endif // OPTIMIZER_CLASS
