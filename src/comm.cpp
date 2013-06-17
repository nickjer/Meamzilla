#include "comm.h"
#include "error.h"
#include <cstdlib>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Comm::Comm() : comm_(MPI_COMM_NULL), root_(0), rank_(-1), size_(0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

Comm::Comm(MPI_Comm comm) : root_(0)
{
  MPI_Comm_dup(comm, &comm_);
  init();
}

/* ---------------------------------------------------------------------- */

Comm::~Comm()
{
  if (*this) MPI_Comm_free(&comm_);
}

/* ----------------------------------------------------------------------
   Initialize rank and size for this comm
------------------------------------------------------------------------- */

void Comm::init()
{
  if (!*this) return;
  MPI_Comm_rank(comm_, &rank_);
  MPI_Comm_size(comm_, &size_);
}

/* ----------------------------------------------------------------------
   return root proc in this communicator
------------------------------------------------------------------------- */

int Comm::get_root() const
{
  return root_;
}

/* ----------------------------------------------------------------------
   return rank of this proc in this communicator
------------------------------------------------------------------------- */

int Comm::get_rank() const
{
  return rank_;
}

/* ----------------------------------------------------------------------
   return number of procs in this communicator
------------------------------------------------------------------------- */

int Comm::get_size() const
{
  return size_;
}

/* ----------------------------------------------------------------------
   determine whether this proc is root or not
------------------------------------------------------------------------- */

bool Comm::is_root() const
{
  return rank_ == root_;
}

/* ----------------------------------------------------------------------
   abort MPI from single proc and exit code
------------------------------------------------------------------------- */

void Comm::abort_one() const
{
  if (*this) MPI_Abort(comm_,1);
  MPI_Finalize();
  std::exit(1);
}

/* ----------------------------------------------------------------------
   abort MPI from all procs and exit code
------------------------------------------------------------------------- */

void Comm::abort_all() const
{
  if (*this) MPI_Barrier(comm_);
  MPI_Finalize();
  std::exit(1);
}

/* ----------------------------------------------------------------------
   MPI split
------------------------------------------------------------------------- */

void Comm::split(int color, int key, Comm *newcomm) const
{
  if (!*this) return;
  MPI_Comm_split( comm_, color, key, &newcomm->comm_ );
  newcomm->init();
}

/* ----------------------------------------------------------------------
   MPI send
------------------------------------------------------------------------- */

int Comm::send(void *buffer, int count, MPI_Datatype datatype, int dest, int tag) const
{
  return MPI_Send( buffer, count, datatype, dest, tag, comm_ );
}

/* ----------------------------------------------------------------------
   MPI send string
------------------------------------------------------------------------- */

int Comm::send(String buffer, int dest, int tag) const
{
  int buff_size = buffer.size();
  send(&buff_size, 1, MPI_INT, dest, tag);
  return send(const_cast<char *>(buffer.c_str()), buff_size, MPI_CHAR, dest, tag);
}

/* ----------------------------------------------------------------------
   MPI recv
------------------------------------------------------------------------- */

int Comm::recv(void *buffer, int count, MPI_Datatype datatype, int src, int tag) const
{
  return MPI_Recv( buffer, count, datatype, src, tag, comm_, MPI_STATUS_IGNORE );
}

/* ----------------------------------------------------------------------
   MPI recv strings
------------------------------------------------------------------------- */

int Comm::recv(String& buffer, int src, int tag) const
{
  int buff_size;
  recv(&buff_size, 1, MPI_INT, src, tag);
  buffer.resize(buff_size);
  return recv(const_cast<char *>(buffer.c_str()), buff_size, MPI_CHAR, src, tag);
}

/* ----------------------------------------------------------------------
   MPI broadcast
------------------------------------------------------------------------- */

void Comm::bcast(void *buffer, int count, MPI_Datatype datatype, int root) const
{
  if (!*this) return;
  MPI_Bcast( buffer, count, datatype, root, comm_ );
}

/* ----------------------------------------------------------------------
   MPI broadcast strings
------------------------------------------------------------------------- */

void Comm::bcast(String& buffer, int root) const
{
  if (!*this) return;

  // Broadcast size of string and resize string
  int buff_size = buffer.size();
  bcast(&buff_size, 1, MPI_INT, root);
  buffer.resize(buff_size);

  // Send string to each proc
  bcast(const_cast<char *>(buffer.c_str()), buff_size, MPI_CHAR, root);
}

/* ----------------------------------------------------------------------
   MPI broadcast list of strings
------------------------------------------------------------------------- */

void Comm::bcast(StringList& buffer, int root) const
{
  if (!*this) return;

  // Broadcast number of lines and resize list of strings for each proc
  int buff_size = buffer.size();
  bcast(&buff_size, 1, MPI_INT, root);
  buffer.resize(buff_size);

  // Send lines to each processor
  for (String& buff_line : buffer)
    bcast(buff_line, root);
}

/* ----------------------------------------------------------------------
   MPI reduce
------------------------------------------------------------------------- */

void Comm::reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root) const
{
  if (!*this) return;
  MPI_Reduce( sendbuf, recvbuf, count, datatype, op, root, comm_ );
}

/* ----------------------------------------------------------------------
   Concatenate strings
------------------------------------------------------------------------- */

void Comm::reduce(String sendbuf, String& recvbuf, int root) const
{
  if (!*this) return;

  recvbuf = sendbuf;

  if (is_root()) {
    for (int i=0; i<size_; ++i) {
      if (i == rank_) continue;
      String tmp_string;
      recv(tmp_string, i);
      recvbuf += tmp_string;
    }
  } else {
    send(sendbuf, root);
  }
}

/* ----------------------------------------------------------------------
   MPI Allreduce
------------------------------------------------------------------------- */

void Comm::all_reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op) const
{
  if (!*this) return;
  MPI_Allreduce( sendbuf, recvbuf, count, datatype, op, comm_ );
}


/* ----------------------------------------------------------------------
   MPI gather
------------------------------------------------------------------------- */

void Comm::gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, int recvcnt, MPI_Datatype recvtype, int root) const
{
  if (!*this) return;
  MPI_Gather( sendbuf, sendcnt, sendtype, recvbuf, recvcnt, recvtype, root, comm_ );
}

/* ----------------------------------------------------------------------
   MPI barrier
------------------------------------------------------------------------- */

void Comm::barrier() const
{
  if (!*this) return;
  MPI_Barrier( comm_ );
}

/* ----------------------------------------------------------------------
   Check if this Comm exists on this processor
------------------------------------------------------------------------- */

Comm::operator bool() const
{
  return (comm_ != MPI_COMM_NULL);
}
