/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifndef MEAMZ_TIMER_H
#define MEAMZ_TIMER_H

namespace MEAMZ_NS
{

class Timer
{
public:
  Timer();
  virtual ~Timer();

  void start();       // Start timer
  void stop();        // Stop timer
  double elapsed();   // Output elapsed time since start
  double get_time();  // Output the value of the time_ variable

protected:

private:
  double time_;

};

} // namespace MEAMZ_NS

#endif // MEAMZ_TIMER_H
