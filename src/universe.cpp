#include "universe.h"
#include "input.h"
#include <fstream>
#include <iostream>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Universe::Universe(class Meamzilla *mmz, MPI_Comm comm) : Pointers(mmz), comm_all(comm), ngroups_(1)
{

}

/* ---------------------------------------------------------------------- */

Universe::~Universe()
{

}

/* ----------------------------------------------------------------------
   set up inter and intra comm's defined by user
------------------------------------------------------------------------- */

void Universe::init()
{
  if (mmz->input)
    mmz->input->parse("ngroups",0,ngroups_);

  if (ngroups_ > comm_all.get_size()) {
    if (comm_all.is_root())
      std::cerr << "ERROR: User defined more groups (ngroups=" << ngroups_ << ") than processors allocated (nprocs=" << comm_all.get_size() << ")" << std::endl;
    abort_all();
  } else if (ngroups_ <= 0) {
    if (comm_all.is_root())
      std::cerr << "ERROR: User defined too few groups (ngroups > 0)" << std::endl;
    abort_all();
  }

  // FIRST DEFINE INTER-COMMUNICATOR

  // Split up all procs into user-defined number of groups
  my_group_ = comm_all.get_rank() / (comm_all.get_size()/ngroups_);
  int key = comm_all.get_rank();
  if (my_group_ >= ngroups_) my_group_ = ngroups_ - 1;  // append remainder to last group
  comm_all.split(my_group_, key, &comm_inter);          // each group will have its own communicator

  // NOW DEFINE INTRA-COMMUNICATOR

  // Get roots from each group and group them together into new comm
  // Ignore the other processors in between
  int color = MPI_UNDEFINED;
  if (comm_inter.is_root()) color = 0;
  comm_all.split(color, key, &comm_intra);


  // COMPUTE NUMBER OF PROCESSORS ON EACH INTER-COMM ON THE GLOBAL ROOT NODE

  if (comm_intra) {
    int nprocs = comm_inter.get_size();

    std::vector<int> nprocs_list (ngroups_, 0);
    comm_intra.gather(&nprocs, 1, MPI_INT, &nprocs_list[0], 1, MPI_INT, comm_intra.get_root());

    if (comm_all.is_root()) {
      int i = 0;
      for (int nprocs : nprocs_list)
        std::cout << "Group #" << i++ << " has " << nprocs << " processors" << std::endl;
    }
  }

  return;
}

/* ----------------------------------------------------------------------
   determine whether current proc is root of universe
------------------------------------------------------------------------- */

bool Universe::is_root() const
{
  return comm_all.is_root();
}

/* ----------------------------------------------------------------------
   return the group number this processor belongs to
------------------------------------------------------------------------- */

int Universe::get_group() const
{
  return my_group_;
}

/* ----------------------------------------------------------------------
   return the total number of groups
------------------------------------------------------------------------- */

int Universe::get_ngroups() const
{
  return ngroups_;
}

/* ----------------------------------------------------------------------
   abort MPI from all procs and exit code
------------------------------------------------------------------------- */

void Universe::abort_all() const
{
  comm_all.abort_all();
  return;
}

/* ----------------------------------------------------------------------
   MPI barrier all procs
------------------------------------------------------------------------- */

void Universe::barrier() const
{
  comm_all.barrier();
  return;
}

/* ----------------------------------------------------------------------
   read a file from main root only and msg contents to all procs
------------------------------------------------------------------------- */

StringList Universe::read_file(const String& filename) const
{
  StringList line_list;

  // Main Root process reads file
  // stores each line in line_list
  if ( comm_all.is_root() ) {
    String line;
    std::ifstream ifs (filename.c_str(), std::ifstream::in);
    if (!ifs) {
      std::cerr << "ERROR: Unable to open file to read: " << filename << std::endl;
      comm_all.abort_one();
    }
    while (!ifs.eof()) {
      std::getline(ifs,line); // get line from file
      if ( line.find("\r") != std::string::npos ) line.erase(line.find("\r"),1);  // remove carriage returns
      for (char character : line) // add to list if find non-whitespace character in line
        if (!std::isspace(character)) { line_list.push_back(line); break; }
    }
    ifs.close();
  }


  // Broadcast line list to all procs
  comm_all.bcast(line_list, comm_all.get_root());

  return line_list;
}
