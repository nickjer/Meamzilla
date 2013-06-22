/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifndef MEAMZ_UNIVERSE_H
#define MEAMZ_UNIVERSE_H

#include "pointers.h"
#include "comm.h"
#include "mpi.h"
#include <string>
#include <vector>

namespace MEAMZ_NS
{

class Universe : protected Pointers
{
public:
  Comm comm_all;
  Comm comm_intra;
  Comm comm_inter;

  Universe(class Meamzilla *, MPI_Comm);
  virtual ~Universe();

  void init();                                  // Initialize inter and intra comm's
  bool is_root() const;                         // Determine whether current proc is root of universe
  int get_group() const;                        // Return the group number this processor belongs to
  int get_ngroups() const;                      // Return the total number of groups
  void abort_all() const;                       // Abort all the cpu's and kill program
  void barrier() const;                         // MPI Barrier all procs
  StringList read_file(const String&) const;    // Read in file and broadcast to all cpu's

protected:

private:
  int ngroups_;                                 // Number of MPI groups (useful for genalg)
  int my_group_;                                // Group number for this processor
};

} // namespace MEAMZ_NS

#endif // MEAMZ_UNIVERSE_H
