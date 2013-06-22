/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifndef MEAMZ_COMM_H
#define MEAMZ_COMM_H

#include "mpi.h"
#include <string>
#include <vector>

namespace MEAMZ_NS
{

typedef std::string String;
typedef std::vector<String> StringList;

class Comm
{
public:
  Comm();
  Comm(MPI_Comm);
  virtual ~Comm();

  int get_root() const;     // Get root proc in this comm
  int get_rank() const;     // Get rank of proc in this comm
  int get_size() const;     // Get # of procs in this comm
  bool is_root() const;     // Whether current proc is root or not
  void abort_one() const;   // Abort MPI from a single processor and exit code
  void abort_all() const;   // Abort MPI from all processors and exit code

  void split(int, int, Comm *) const;                                             // MPI Split
  int send(void *, int, MPI_Datatype, int, int = 0) const;                        // MPI Send
  int send(String, int, int = 0) const;                                           // MPI Send string
  int recv(void *, int, MPI_Datatype, int, int = 0) const;                        // MPI Recv
  int recv(String&, int, int = 0) const;                                          // MPI Recv string
  void bcast(void *, int, MPI_Datatype, int) const;                               // MPI Broadcast
  void bcast(String&, int) const;                                                 // MPI Broadcast strings
  void bcast(StringList&, int) const;                                             // MPI Broadcast list of strings
  void reduce(void *, void *, int, MPI_Datatype, MPI_Op, int) const;              // MPI Reduce
  void reduce(String, String&, int) const;                                        // Concatenate strings
  void all_reduce(void *, void *, int, MPI_Datatype, MPI_Op) const;               // MPI Allreduce
  void gather(void *, int, MPI_Datatype, void *, int, MPI_Datatype, int) const;   // MPI gather
  void barrier() const;                                                           // MPI Barrier

  explicit operator bool() const;    // Check if this Comm exists on this processor

protected:
private:
  MPI_Comm comm_;
  int root_;
  int rank_;
  int size_;

  void init();              // Initialize rank and size for this comm
};

} // namespace MEAMZ_NS

#endif // MEAMZ_COMM_H
