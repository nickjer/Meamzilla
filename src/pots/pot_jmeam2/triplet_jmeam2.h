#ifdef TRIPLET_CLASS

TripletStyle(jmeam2,TripletJMEAM2)

#else

#ifndef MEAMZ_TRIPLET_JMEAM2_H
#define MEAMZ_TRIPLET_JMEAM2_H

#include "../pot_jmeam/triplet_jmeam.h"

namespace MEAMZ_NS
{

class TripletJMEAM2 : public TripletJMEAM
{
public:
  int q2_idx;                        // Index for q- alloy potential used for this pair
  int q2_knot;                       // -1 if outside range of q- potential, anything else otherwise

  double q2;                         // q-potential value for triplet
  double dq2;                        // Derivative of q- potential for this triplet: d(q)/dcos

  TripletJMEAM2();
  virtual ~TripletJMEAM2();

  virtual TripletJMEAM2* clone() const;              // Copy polymorphic objects

  virtual bool check_triplet(Atom *, Potential *);  // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_TRIPLET_JMEAM2_H
#endif // TRIPLET_CLASS
