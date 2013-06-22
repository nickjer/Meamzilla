/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(jmeam2,AtomJMEAM2)

#else

#ifndef MEAMZ_ATOM_JMEAM2_H
#define MEAMZ_ATOM_JMEAM2_H

#include "../pot_jmeam/atom_jmeam.h"

namespace MEAMZ_NS
{

class AtomJMEAM2 : public AtomJMEAM
{
public:
  AtomJMEAM2();
  virtual ~AtomJMEAM2();

  virtual AtomJMEAM2* clone() const;  // Copy polymorphic objects

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_JMEAM2_H
#endif // ATOM_CLASS
