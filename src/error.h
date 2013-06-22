/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifndef MEAMZ_ERROR_H
#define MEAMZ_ERROR_H

#include <string>
#include <stdexcept>

namespace MEAMZ_NS
{

class Error : public std::runtime_error
{
public:
  Error(const std::string& msg) : std::runtime_error(msg) {}
};

} // namespace MEAMZ_NS

#endif // MEAMZ_ERROR_H
