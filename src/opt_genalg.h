/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef OPTIMIZER_CLASS

OptimizerStyle(genalg,OptGenAlg)

#else

#ifndef MEAMZ_OPT_GENALG_H
#define MEAMZ_OPT_GENALG_H

#include "optimizer.h"

namespace MEAMZ_NS
{

class OptGenAlg : public Optimizer
{
public:
  OptGenAlg(class Meamzilla *, std::ostream * = nullptr);
  virtual ~OptGenAlg();

  virtual void init();              // Initialize the optimizer
  virtual void compute(PotList&);   // Use the optimizer on a potential list

protected:
  int pop_size_;          // Population size
  double cross_rate_;     // Crossover rate
  double mut_rate_;       // Mutation rate
  double init_scale_;     // Initial fluctutation in scale of height of potential parameters in pop.
  double fit_rate_;       // Rate of how long to breed fittest potentials before moving to next fittest pot
  double rescale_rate_;   // Rate at which we rescale the potentials at each genalg step
  int order_breed_;       // Do we breed best potentials in order or breed them randomly?
  int gen_save_;          // Number of top potentials to save

  int num_powell_;        // Number of Powell minimization steps

  // Constants used during simulated annealing
  static constexpr double EPS = 0.001;
  static constexpr int NEPS = 10;

  PotList compute_potlist(const PotList&, int) const;

private:

};

} // MEAMZ_NS

#endif // MEAMZ_OPT_GENALG_H
#endif // OPTIMIZER_CLASS
