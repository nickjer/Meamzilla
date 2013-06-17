#ifndef MEAMZ_INPUT_H
#define MEAMZ_INPUT_H

#include "pointers.h"
#include "universe.h"
#include <sstream>
#include <iostream>

namespace MEAMZ_NS
{

class Input : protected Pointers
{
public:
  Input(class Meamzilla *, const StringList);
  virtual ~Input();

  template <typename... Args> void parse(const String, int, Args&...);

  void parse_line(std::istringstream& line_str) {}        // termination of recursive parse_line fn

  template <typename Arg1, typename... Args>
  void parse_line(std::istringstream&, Arg1&, Args&...) throw(int);

  void parse_stringlist(const String, int, StringList&);

protected:
private:
    StringList input_params_;         // User input from parameters file
};

/* ----------------------------------------------------------------------
   Find the parameter in parameter list and parse the line
------------------------------------------------------------------------- */

template <typename... Args>
void Input::parse(const String given_param, int must_be_found, Args&... args)
{
  int success = 0; // Marker whether we found the given param

  for (String line : input_params_) {
    std::istringstream line_str (line, std::istringstream::in);

    String read_param;
    line_str >> read_param; // Read in parameter from file

    if (read_param == given_param) {  // Found
      try {
        parse_line(line_str, args...);
      } catch(int error) {
        if (mmz->universe->comm_all.is_root())
          std::cerr << "ERROR: Parameter line below is ill-defined:\n" << line << std::endl;
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

/* ----------------------------------------------------------------------
   Recursively parse each argument
------------------------------------------------------------------------- */

template <typename Arg1, typename... Args>
void Input::parse_line(std::istringstream& line_str, Arg1& arg1, Args&... args) throw(int)
{
  line_str >> arg1; // Read in argument
  if (line_str.fail()) throw 1; // Throw error if stream fails
  parse_line(line_str, args...);  // Recursively do the same for the remaining arguments
  return;
}

} // namespace MEAMZ_NS

#endif // MEAMZ_INPUT_H
