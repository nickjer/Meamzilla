#ifdef POTENTIAL_CLASS

PotentialStyle(pair/spline,PotPairSpline)

#else

#ifndef MEAMZ_POT_PAIR_SPLINE_H
#define MEAMZ_POT_PAIR_SPLINE_H

#include "../pot_pair/pot_pair.h"

namespace MEAMZ_NS {

class PotPairSpline : public PotPair
{
public:
  PotPairSpline(class Meamzilla *, int);
  virtual ~PotPairSpline();

  virtual PotPairSpline* clone() const;                       // Copy polymorphic objects

  virtual double compute(const Comm&, ErrorVec * = nullptr);    // Compute error sum/vector

protected:

private:
  template <class T>
  double fast_compute(const Comm&, ErrorVec *);

};

} // namespace MEAMZ_NS

#endif // MEAMZ_POT_PAIR_SPLINE_H
#endif // POTENTIAL_CLASS
