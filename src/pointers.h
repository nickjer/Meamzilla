/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifndef MEAMZ_POINTERS_H
#define MEAMZ_POINTERS_H

#include "meamzilla.h"

namespace MEAMZ_NS
{

class Pointers
{
public:
  Pointers(Meamzilla *ptr) : mmz(ptr) {}
  virtual ~Pointers() {}

protected:
  Meamzilla *mmz;

private:
};

} // namespace MEAMZ_NS

#endif // MEAMZ_POINTERS_H
