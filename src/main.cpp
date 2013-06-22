/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "meamzilla.h"
#include "mpi.h"


using namespace MEAMZ_NS;

int main(int narg, char **arg)
{
  MPI_Init(&narg, &arg);

  Meamzilla *mmz = new Meamzilla(narg, arg, MPI_COMM_WORLD);

  mmz->run();

  delete mmz;

  MPI_Finalize();

  return 0;
}
