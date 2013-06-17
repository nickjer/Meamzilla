#ifdef TRIPLET_CLASS

TripletStyle(jmeam2/spline,TripletJMEAM2Spline)

#else

#ifndef MEAMZ_TRIPLET_JMEAM2_SPLINE_H
#define MEAMZ_TRIPLET_JMEAM2_SPLINE_H

#include "../pot_jmeam2/triplet_jmeam2.h"

namespace MEAMZ_NS
{

class TripletJMEAM2Spline : public TripletJMEAM2
{
public:
  double g_shift;
  double q_shift;
  double q2_shift;

  TripletJMEAM2Spline();
  virtual ~TripletJMEAM2Spline();

  virtual TripletJMEAM2Spline* clone() const;        // Copy polymorphic objects

  virtual bool check_triplet(Atom *, Potential *);  // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_TRIPLET_JMEAM2_SPLINE_H
#endif // TRIPLET_CLASS
