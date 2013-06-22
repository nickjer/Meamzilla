/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef TRIPLET_CLASS

TripletStyle(jmeam,TripletJMEAM)

#else

#ifndef MEAMZ_TRIPLET_JMEAM_H
#define MEAMZ_TRIPLET_JMEAM_H

#include "../pot_meam/triplet_meam.h"

namespace MEAMZ_NS
{

class TripletJMEAM : public TripletMEAM
{
public:
  int q_idx;                        // Index for q- alloy potential used for this pair
  int q_knot;                       // -1 if outside range of q- potential, anything else otherwise

  double q;                         // q-potential value for triplet
  double dq;                        // Derivative of q- potential for this triplet: d(q)/dcos

  TripletJMEAM();
  virtual ~TripletJMEAM();

  virtual TripletJMEAM* clone() const;              // Copy polymorphic objects

  virtual bool check_triplet(Atom *, Potential *);  // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_TRIPLET_JMEAM_H
#endif // TRIPLET_CLASS
