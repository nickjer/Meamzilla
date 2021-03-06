/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(meam,AtomMEAM)

#else

#ifndef MEAMZ_ATOM_MEAM_H
#define MEAMZ_ATOM_MEAM_H

#include "../pot_eam/atom_eam.h"

namespace MEAMZ_NS
{

class AtomMEAM : public AtomEAM
{
public:
  AtomMEAM();
  virtual ~AtomMEAM();

  virtual AtomMEAM* clone() const;  // Copy polymorphic objects

  virtual int cmp_value() const;    // Value used when sorting atoms

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_MEAM_H
#endif // ATOM_CLASS
