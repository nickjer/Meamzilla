#include "input.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Input::Input(class Meamzilla *mmz, const StringList input_params) : Pointers(mmz), input_params_(input_params)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

Input::~Input()
{
  //dtor
}

/* ----------------------------------------------------------------------
   Find the parameter in parameter list and parse the line
------------------------------------------------------------------------- */

void Input::parse_stringlist(const String given_param, int must_be_found, StringList& args)
{
  int success = 0; // Marker whether we found the given param

  for (String line : input_params_) {
    std::istringstream line_str (line, std::istringstream::in);

    String read_param;
    line_str >> read_param; // Read in parameter from file

    if (read_param == given_param) {  // Found
      while (!line_str.eof()) {
        String tmpstring;
        line_str >> tmpstring;
        if (line_str.fail()) break;
        args.push_back(tmpstring);
      }
      if (!args.size() && must_be_found) {
        if (mmz->universe->comm_all.is_root())
          std::cerr << "ERROR: Parameter \"" << given_param << "\" is ill-defined" << std::endl;
        mmz->universe->abort_all();
      }
      success = 1;  // Mark as found and leave loop
      break;
    }
  }

  if (!success && must_be_found) {
    if (mmz->universe->comm_all.is_root())
      std::cerr << "ERROR: Parameter not found in file: \"" << given_param << "\"" << std::endl;
    mmz->universe->abort_all();
  }
  return;
}
