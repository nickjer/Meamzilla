/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifndef MEAMZ_OPT_LIST_H
#define MEAMZ_OPT_LIST_H

#include "pointers.h"
#include "optimizer.h"

namespace MEAMZ_NS
{

class OptList : protected Pointers
{
public:
  OptList(class Meamzilla *);
  virtual ~OptList();

  void init();                  // Initialize this class
  void compute(PotList&);       // Compute list of optimizers using list of potentials

protected:

private:
  int nopts_;                       // Number of potentials
  std::vector<Optimizer *> opts_;   // List of optimizers

};

} // namespace MEAMZ_NS

#endif // MEAMZ_OPT_LIST_H
