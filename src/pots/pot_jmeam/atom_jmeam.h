/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(jmeam,AtomJMEAM)

#else

#ifndef MEAMZ_ATOM_JMEAM_H
#define MEAMZ_ATOM_JMEAM_H

#include "../pot_meam/atom_meam.h"

namespace MEAMZ_NS
{

class AtomJMEAM : public AtomMEAM
{
public:
  AtomJMEAM();
  virtual ~AtomJMEAM();

  virtual AtomJMEAM* clone() const;  // Copy polymorphic objects

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_JMEAM_H
#endif // ATOM_CLASS
