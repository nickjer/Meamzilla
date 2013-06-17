#ifdef TRIPLET_CLASS

TripletStyle(jmeam/spline,TripletJMEAMSpline)

#else

#ifndef MEAMZ_TRIPLET_JMEAM_SPLINE_H
#define MEAMZ_TRIPLET_JMEAM_SPLINE_H

#include "../pot_jmeam/triplet_jmeam.h"

namespace MEAMZ_NS
{

class TripletJMEAMSpline : public TripletJMEAM
{
public:
  double g_shift;
  double q_shift;

  TripletJMEAMSpline();
  virtual ~TripletJMEAMSpline();

  virtual TripletJMEAMSpline* clone() const;        // Copy polymorphic objects

  virtual bool check_triplet(Atom *, Potential *);  // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_TRIPLET_JMEAM_SPLINE_H
#endif // TRIPLET_CLASS
