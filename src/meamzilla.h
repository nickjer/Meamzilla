#ifndef MEAMZ_MEAMZILLA_H
#define MEAMZ_MEAMZILLA_H

#include "mpi.h"
#include <vector>
#include <string>

namespace MEAMZ_NS
{

typedef std::string String;
typedef std::vector<String> StringList;

class Meamzilla
{
public:
  class Universe *universe;        // Universe of processors
  class Input *input;              // Input defined by user (usually in a parameter file)
  class Timer *timer;              // Timer used for benchmarking
  class Random *random;            // Random number generator
  class Config *config;            // Configurations of supercells
  class AtomVec *atomvec;          // List of atoms for force computations (different for each processor)

  class PotList *potlist;          // List of potentials that need to be optimized
  class OptList *optlist;          // List of optimizers to be run in sequential order

  Meamzilla(MPI_Comm);
  Meamzilla(int, char **, MPI_Comm);
  virtual ~Meamzilla();

  void init();      // Initialize members of this class
  void run();       // Run Meamzilla

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_MEAMZILLA_H
