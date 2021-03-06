/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef POTENTIAL_CLASS

PotentialStyle(jmeam/spline,PotJMEAMSpline)

#else

#ifndef MEAMZ_POT_JMEAM_SPLINE_H
#define MEAMZ_POT_JMEAM_SPLINE_H

#include "../pot_jmeam/pot_jmeam.h"

namespace MEAMZ_NS {

class PotJMEAMSpline : public PotJMEAM
{
public:
  PotJMEAMSpline(class Meamzilla *, int);
  virtual ~PotJMEAMSpline();

  virtual PotJMEAMSpline* clone() const;                     // Copy polymorphic objects

  virtual double compute(const Comm&, ErrorVec * = nullptr);      // Compute error sum/vector

protected:
  virtual void compute_densities(const Comm&, Vector<double> *);  // Compute densities for each atom

private:
  template <class T>
  double fast_compute(const Comm&, ErrorVec *);

  template <class T>
  void fast_compute_densities(const Comm&, Vector<double> *);

};

} // namespace MEAMZ_NS

#endif // MEAMZ_POT_JMEAM_SPLINE_H
#endif // POTENTIAL_CLASS
