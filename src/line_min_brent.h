/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifndef MEAMZ_LINE_MIN_BRENT_H
#define MEAMZ_LINE_MIN_BRENT_H

#include "line_min.h"

namespace MEAMZ_NS
{

class LineMinBrent : public LineMin
{
public:
  LineMinBrent();
  virtual ~LineMinBrent();

  double find_minima();   // Find the minima along line, return error_sum

protected:
  static constexpr unsigned int ITMAX = 100;      // Max # allowed iterations
  static constexpr double TOL = 1.0e-1;           // Tolerance to know when to stop
  static constexpr double ZEPS = 1.0e-9;          // Small # that protects against trying to achieve fractional accuracy
                                                  // for a minimum that happens to be exactly 0
private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_LINE_MIN_BRENT_H
