/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS

AtomStyle(eam,AtomEAM)

#else

#ifndef MEAMZ_ATOM_EAM_H
#define MEAMZ_ATOM_EAM_H

#include "../pot_pair/atom_pair.h"

namespace MEAMZ_NS
{

class AtomEAM : public AtomPair
{
public:
  int F_idx;              // index of embedding function for this atom (equal to atom type)

  AtomEAM();
  virtual ~AtomEAM();

  virtual AtomEAM* clone() const;  // Copy polymorphic objects

  virtual void setup_atom_pot();   // Setup potential specific atom properties

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_EAM_H
#endif // ATOM_CLASS
