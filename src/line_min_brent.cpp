#include "line_min_brent.h"
#include "error.h"

#include <cmath>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

LineMinBrent::LineMinBrent() : LineMin()
{
  //ctor
}

/* ---------------------------------------------------------------------- */

LineMinBrent::~LineMinBrent()
{
  //dtor
}

/* ----------------------------------------------------------------------
   find the minima along line, return error_sum
------------------------------------------------------------------------- */

double LineMinBrent::find_minima()
{
  double x_left = (xa_ < xc_ ? xa_ : xc_);    // Left bound
  double x_right = (xa_ > xc_ ? xa_ : xc_);   // Right bound
  double v = xb_;                             // Previous value of w
  double w = v;                               // Point with 2nd least fn value
  double x = v;                               // Point with least fn value
  double d = 0.;                              // Distance moved on current step
  double e = 0.;                              // This will be the distance moved on the step before last

  double Fx = Fb_;
  double Fv = Fx;
  double Fw = Fx;

  for (unsigned int count=0; count<ITMAX; ++count) {
    double u = 0.;  // Temporary x-position that may lead to lowest fn value
    double xm = 0.5*(x_left + x_right);   // Get midpoint
    double tol1 = TOL * std::abs(x) + ZEPS;
    double tol2 = 2. * tol1;

    // SUCCESS: we need to close up here
    if (std::abs(x - xm) <= (tol2 - 0.5*(x_right - x_left))) {
      xmin1_ = x; // min
      xmin2_ = w; // min2
      Fxmin1_ = Fx;
      Fxmin2_ = Fw;
      return Fxmin1_;
    }

    if (std::abs(e) > tol1) {
      // Construct a trial parabolic fit.
      double r = (x-w)*(Fx-Fv);
      double q = (x-v)*(Fx-Fw);
      double p = (x-v)*q - (x-w)*r;
      q = 2.*(q-r);
      if (q>0.) p = -p;
      else q = -q;
      double etemp = e;
      e = d;
      // Determine acceptability of parabolic fit
      if (std::abs(p) >= std::abs(0.5*q*etemp) || p <= q*(x_left-x) || p >= q*(x_right-x)) {
        // Fail: Take a golden section step into the larger of the two segments
        e = (x < xm) ? x_right - x : -(x - x_left);
        d = CGOLD * e;  // Golden section step
      } else {
        // Parabolic interpolation is attempted, fitting through the points x, v, and w.
        d = p/q;  // Take the parabolic step
        u = x+d;
        // Make sure step isn't outside bounds of bracket
        if ((u - x_left) < tol2 || (x_right - u) < tol2)
          d = (x < xm) ? tol1 : -tol1;  // Determine new direction
      }
    } else {
      // Take golden section step
      e = (x < xm) ? x_right - x : -(x - x_left);
      d = CGOLD * e;  // Golden section step
    }
    // Arrive here with d computed either from parabolic fit, or else from golden section.
    if (std::abs(d) >= tol1) u = x+d;
    else u = x + ((d > 0) ? tol1 : -tol1);
    // This is the one function evaluation per iteration,and now we have to decide what to do with our function evaluation
    double Fu = eval_(u);
    if (Fu <= Fx) {
      if (u >= x) x_left = x;
      else x_right = x;
      v = w;
      Fv = Fw;
      w = x;
      Fw = Fx;
      x = u;
      Fx = Fu;
    } else {
      if (u < x) x_left = u;
      else x_right = u;
      if (Fu <= Fw || w == x) {
        v = w;
        Fv = Fw;
        w = u;
        Fw = Fu;
      } else if (Fu <= Fv || v == x || v == w) {
        v = u;
        Fv = Fu;
      }
    }
  }

  throw Error("ERROR: too many iterations in Brent minimization - ");

  return 0.;
}
