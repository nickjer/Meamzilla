/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "meamzilla.h"
#include "universe.h"
#include "input.h"
#include "timer.h"
#include "random.h"
#include "config.h"
#include "atom_vec.h"
#include "pot_list.h"
#include "opt_list.h"
#include <iostream>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Meamzilla::Meamzilla(MPI_Comm communicator)
{
  universe = new Universe(this, communicator);
  timer = new Timer();
  random = new Random(this);
  config = new Config(this);
  atomvec = new AtomVec(this);

  potlist = new PotList(this);
  optlist = new OptList(this);
}

/* ---------------------------------------------------------------------- */

Meamzilla::Meamzilla(int narg, char **arg, MPI_Comm communicator) : Meamzilla(communicator)
{
  String paramsfilename;

  String arg_error = String("ERROR: Invalid command-line argument\n") + String(arg[0]) + String(" -p <parameter file>");

  // Read in command line arguments
  int iarg = 1;
  while (iarg < narg) {
    String command = arg[iarg];

    if (command == "-param" || command == "-p") {
      if (iarg+2 > narg) {
        if (universe->comm_all.is_root())
          std::cerr << arg_error << std::endl;
        universe->abort_all();
      }
      paramsfilename = arg[iarg + 1];
      iarg += 2;
    } else {
      if (universe->comm_all.is_root())
        std::cerr << arg_error << std::endl;
      universe->abort_all();
    }
  }

  // Need parameter file
  if (paramsfilename == "") {
    if (universe->comm_all.is_root())
      std::cerr << arg_error << std::endl;
    universe->abort_all();
  }

  // Read in parameter file and store it in input
  input = new Input(this, universe->read_file(paramsfilename));

  init();
}

/* ---------------------------------------------------------------------- */

Meamzilla::~Meamzilla()
{
  if (input) delete input;
  delete atomvec;
  delete config;
  delete potlist;
  delete random;
  delete universe;
}

/* ----------------------------------------------------------------------
   initialize members of this class
------------------------------------------------------------------------- */

void Meamzilla::init()
{
  universe->init();   // initialize universe as defined by user (break up into groups)
  random->init();     // initialize random number generators
  potlist->init();    // initialize list of potentials
  optlist->init();    // initialize list of optimizers
  config->init();     // initialize configurations (must be done after potlist to get atom types)
  atomvec->init();    // initialize list of atoms to be computed for each processor (must be done after config)
}

/* ----------------------------------------------------------------------
   run Meamzilla
------------------------------------------------------------------------- */

void Meamzilla::run()
{
  PotList& plist = *potlist;  // for quick reference

  // Output total sum error
  {
    std::ostringstream oss;
    for (int i=0; i<plist.get_group_npots(); ++i) {
      double F = plist[i]->compute(universe->comm_inter);
      if (universe->comm_inter.is_root()) {
        oss << "G#" << universe->get_group() << " ";
        oss << "P#" << plist[i]->get_global_idx() << " - ";
        oss << "Initial total sum error = ";
        oss << std::fixed << F << std::endl;
      }
    }
    std::string final_string;
    universe->comm_intra.reduce(oss.str(), final_string, universe->comm_intra.get_root());
    if (universe->comm_intra.is_root())
      std::cout << final_string << std::endl;
    universe->barrier();
  }

  // Rescale potential in beginning if necessary
  if (universe->is_root()) {
    if (potlist->get_rescale() == 2)
      std::cout << "Checking if rescaling is necessary..." << std::endl;
    else if (potlist->get_rescale() == 3)
      std::cout << "Rescaling potential..." << std::endl;
  }
  universe->barrier();

  for (int i=0; i<plist.get_group_npots(); ++i)
    plist[i]->rescale(universe->comm_inter, &std::cout, potlist->get_rescale() - 2);
  universe->barrier();

  optlist->compute(*potlist);

  plist.write_pots(plist.get_end_file());
  plist.write_lmps(plist.get_lmp_file());
  plist.write_data(plist.get_dat_file());

  // Output total sum error
  {
    std::ostringstream oss;
    for (int i=0; i<plist.get_group_npots(); ++i) {
      double F = plist[i]->compute(universe->comm_inter);
      if (universe->comm_inter.is_root()) {
        oss << "G#" << universe->get_group() << " ";
        oss << "P#" << plist[i]->get_global_idx() << " - ";
        oss << "Final total sum error = ";
        oss << std::fixed << F << std::endl;
      }
    }
    std::string final_string;
    universe->comm_intra.reduce(oss.str(), final_string, universe->comm_intra.get_root());
    if (universe->comm_intra.is_root())
      std::cout << std::endl << final_string << std::endl;
    universe->barrier();
  }

  return;
}
