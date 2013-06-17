#ifdef PAIR_CLASS

PairStyle(jmeam2/spline,PairJMEAM2Spline)

#else

#ifndef MEAMZ_PAIR_JMEAM2_SPLINE_H
#define MEAMZ_PAIR_JMEAM2_SPLINE_H

#include "../pot_jmeam2/pair_jmeam2.h"

namespace MEAMZ_NS
{

class PairJMEAM2Spline : public PairJMEAM2
{
public:
  double phi_shift;
  double rho_shift;
  double f_shift;
  double p_shift;
  double p2_shift;

  PairJMEAM2Spline();
  virtual ~PairJMEAM2Spline();

  virtual PairJMEAM2Spline* clone() const;      // Copy polymorphic objects

  virtual bool check_pair(Atom *, Potential *); // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_PAIR_JMEAM2_H
#endif // PAIR_CLASS
