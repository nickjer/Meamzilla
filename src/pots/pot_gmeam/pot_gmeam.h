/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef POTENTIAL_CLASS

PotentialStyle(gmeam,PotGMEAM)

#else

#ifndef MEAMZ_POT_GMEAM_H
#define MEAMZ_POT_GMEAM_H

#include "../pot_meam/pot_meam.h"

namespace MEAMZ_NS {

class PotGMEAM : public PotMEAM
{
public:
  PotGMEAM(class Meamzilla *, int);
  virtual ~PotGMEAM();

  virtual PotGMEAM* clone() const;                    // Copy polymorphic objects

  virtual double compute(const Comm&, ErrorVec * = nullptr);      // Compute error sum/vector

protected:
  virtual void write_punish(const ErrorVec&, String, int) const;  // Write punishment/constraint errors (none for this pot type)

  virtual int rescale_3body(const Comm&, std::ostream *, int);    // Rescale MEAM part of the potential: f*f*g

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_POT_GMEAM_H
#endif // POTENTIAL_CLASS
