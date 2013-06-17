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
