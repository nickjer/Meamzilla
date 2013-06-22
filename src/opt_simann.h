/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef OPTIMIZER_CLASS

OptimizerStyle(simann,OptSimAnn)

#else

#ifndef MEAMZ_OPT_SIMANN_H
#define MEAMZ_OPT_SIMANN_H

#include "optimizer.h"

namespace MEAMZ_NS
{

class OptSimAnn : public Optimizer
{
public:
  OptSimAnn(class Meamzilla *, std::ostream * = nullptr);
  virtual ~OptSimAnn();

  virtual void init();              // Initialize the optimizer
  virtual void compute(PotList&);   // Use the optimizer on a potential list

  virtual Potential* compute_pot(const Potential *, int) const;  // Optimization loop

protected:
  double anneal_temp_;              // Annealing temperature

  // Constants used during simulated annealing
  static constexpr double EPS = 0.1;
  static constexpr int NEPS = 4;
  //static const int NSTEP = 2;
  static constexpr int NSTEP = 20;
  static constexpr double STEPVAR = 2.0;
  static constexpr double TEMPVAR = 0.85;

  void print_simann(int, double, int, double, double, int, int) const;      // Print results
  int metropolis(Potential*&, Potential *, double *, double,                // Do a metropolis step (return 1 if accepted)
                 double, Potential*&, double *) const;

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_OPT_SIMANN_H
#endif // OPTIMIZER_CLASS
