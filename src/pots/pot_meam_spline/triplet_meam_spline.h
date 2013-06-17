#ifdef TRIPLET_CLASS

TripletStyle(meam/spline,TripletMEAMSpline)

#else

#ifndef MEAMZ_TRIPLET_MEAM_SPLINE_H
#define MEAMZ_TRIPLET_MEAM_SPLINE_H

#include "../pot_meam/triplet_meam.h"

namespace MEAMZ_NS
{

class TripletMEAMSpline : public TripletMEAM
{
public:
  double g_shift;

  TripletMEAMSpline();
  virtual ~TripletMEAMSpline();

  virtual TripletMEAMSpline* clone() const;               // Copy polymorphic objects

  virtual bool check_triplet(Atom *, Potential *);  // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_TRIPLET_MEAM_SPLINE_H
#endif // TRIPLET_CLASS
