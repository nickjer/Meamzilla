/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef POTENTIAL_CLASS

PotentialStyle(jmeam,PotJMEAM)

#else

#ifndef MEAMZ_POT_JMEAM_H
#define MEAMZ_POT_JMEAM_H

#include "../pot_meam/pot_meam.h"

namespace MEAMZ_NS {

class PotJMEAM : public PotMEAM
{
public:
  PotJMEAM(class Meamzilla *, int);
  virtual ~PotJMEAM();

  virtual PotJMEAM* clone() const;                    // Copy polymorphic objects

  virtual void check_if_similar(Potential *) const;   // Check that a potential is similar enough to this
                                                      // potential for pairs/triplets to have identical info

  virtual double compute(const Comm&, ErrorVec * = nullptr);      // Compute error sum/vector

protected:
  virtual void write_punish(const ErrorVec&, String, int) const;  // Write punishment/constraint errors (none for this pot type)

  virtual int rescale_3body(const Comm&, std::ostream *, int);    // Rescale MEAM part of the potential: f*f*g

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_POT_JMEAM_H
#endif // POTENTIAL_CLASS
