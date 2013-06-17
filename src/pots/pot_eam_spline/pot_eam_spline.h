#ifdef POTENTIAL_CLASS

PotentialStyle(eam/spline,PotEAMSpline)

#else

#ifndef MEAMZ_POT_EAM_SPLINE_H
#define MEAMZ_POT_EAM_SPLINE_H

#include "../pot_eam/pot_eam.h"

namespace MEAMZ_NS {

class PotEAMSpline : public PotEAM
{
public:
  PotEAMSpline(class Meamzilla *, int);
  virtual ~PotEAMSpline();

  virtual PotEAMSpline* clone() const;                      // Copy polymorphic objects

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

#endif // MEAMZ_POT_EAM_SPLINE_H
#endif // POTENTIAL_CLASS
