#ifdef POTENTIAL_CLASS

PotentialStyle(jmeam2/spline,PotJMEAM2Spline)

#else

#ifndef MEAMZ_POT_JMEAM2_SPLINE_H
#define MEAMZ_POT_JMEAM2_SPLINE_H

#include "../pot_jmeam2/pot_jmeam2.h"

namespace MEAMZ_NS {

class PotJMEAM2Spline : public PotJMEAM2
{
public:
  PotJMEAM2Spline(class Meamzilla *, int);
  virtual ~PotJMEAM2Spline();

  virtual PotJMEAM2Spline* clone() const;                     // Copy polymorphic objects

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

#endif // MEAMZ_POT_JMEAM2_SPLINE_H
#endif // POTENTIAL_CLASS
