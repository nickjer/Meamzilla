/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "line_min.h"
#include "error.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

LineMin::LineMin()
  : xa_(0.0), xb_(0.0), xc_(0.0), Fa_(0.0), Fb_(0.0), Fc_(0.0),
    xmin1_(0.0), xmin2_(0.0), Fxmin1_(0.0), Fxmin2_(0.0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

LineMin::~LineMin()
{
  //dtor
}

/* ----------------------------------------------------------------------
   get xmin1
------------------------------------------------------------------------- */

double LineMin::get_xmin1() const
{
  return xmin1_;
}

/* ----------------------------------------------------------------------
   get xmin2
------------------------------------------------------------------------- */

double LineMin::get_xmin2() const
{
  return xmin2_;
}

/* ----------------------------------------------------------------------
   get Fmin1
------------------------------------------------------------------------- */

double LineMin::get_Fmin1() const
{
  return Fxmin1_;
}

/* ----------------------------------------------------------------------
   get Fmin2
------------------------------------------------------------------------- */

double LineMin::get_Fmin2() const
{
  return Fxmin2_;
}

/* ----------------------------------------------------------------------
   set function to be minimized
------------------------------------------------------------------------- */

void LineMin::set_func(std::function<double(double)>& eval)
{
  eval_ = eval;
  get_bracket();
}

/* ----------------------------------------------------------------------
   bracket the minimum (xa,xb,xc) where f(xa)>f(xb)<f(xc)
------------------------------------------------------------------------- */

void LineMin::get_bracket()
{
  double x_left = 0.0;
  double x_middle;
  double x_right = 0.1;
  double F_left = eval_(x_left);
  double F_middle;
  double F_right = eval_(x_right);  // evaluate at x_right

  // Define the middle point now using golden sections
  if (F_right >= F_left) {  // most likely min is to the left
    x_middle = x_left;
    F_middle = F_left;
    x_left = -(x_right - x_middle)/CGOLD + x_right; // left side is approx 62% and right side is 38% of total new length
    F_left = eval_(x_left);
  } else {  // most likely min is to the right
    x_middle = x_right;
    F_middle = F_right;
    x_right = (x_middle - x_left)/CGOLD + x_left; // left side is 38% and right side is 62% of total new length
    F_right = eval_(x_right);
  }

  for (unsigned int count=0; count<MAX_BRACKET; ++count) {
    // We need to check that F_left > F_middle < F_right
    if (F_middle < F_left) {
      if (F_middle < F_right) {
        // SUCCESS!!! The minima is inside of our bracket
        xa_ = x_left;
        xb_ = x_middle;
        xc_ = x_right;
        Fa_ = F_left;
        Fb_ = F_middle;
        Fc_ = F_right;
        return;
      } else if (F_middle > F_right) {
        // Just go right
        x_left = x_middle;
        F_left = F_middle;
        x_middle = x_right;
        F_middle = F_right;
        x_right = (x_middle - x_left)/CGOLD + x_left; // left side is 38% and right side is 62% of total new length
        F_right = eval_(x_right);
      } else {  // F_center == F_right
        // Pathological: Search between middle and right
        // This means a change from original algorithm
        // But I doubt this ever occurs
        x_right = (x_right - x_left) * CGOLD + x_right; // the three points are changed a bit from orig. percentage
        F_right = eval_(x_right);
      }
    } else if (F_middle > F_left) {
      // Just go left
      // Just go to the left
      x_right = x_middle;
      F_right = F_middle;
      x_middle = x_left;
      F_middle = F_left;
      x_left = -(x_right - x_middle)/CGOLD + x_right; // left side is approx 62% and right side is 38% of total new length
      F_left = eval_(x_left);
    } else {  // F_center == F_left
      // Pathological: Search between middle and right
      // This means a change from original algorithm
      // But I doubt this ever occurs
      x_left = -(x_right - x_left) * CGOLD + x_left; // the three points are changed a bit from orig. percentage
      F_left = eval_(x_left);
    }
  }

  throw Error("ERROR: problems with bracketing of minimum - ");

  return;
}
