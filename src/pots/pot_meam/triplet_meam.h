#ifdef TRIPLET_CLASS

TripletStyle(meam,TripletMEAM)

#else

#ifndef MEAMZ_TRIPLET_MEAM_H
#define MEAMZ_TRIPLET_MEAM_H

#include "../../triplet.h"

namespace MEAMZ_NS
{

class TripletMEAM : public Triplet
{
public:
  int g_idx;                        // Index for g- alloy potential used for this pair
  int g_knot;                       // -1 if outside range of g- potential, anything else otherwise

  double g;                         // g-potential value for triplet
  double dg;                        // Derivative of g- potential for this triplet: d(g)/dcos

  TripletMEAM();
  virtual ~TripletMEAM();

  virtual TripletMEAM* clone() const;               // Copy polymorphic objects

  virtual bool check_triplet(Atom *, Potential *);  // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_TRIPLET_MEAM_H
#endif // TRIPLET_CLASS
