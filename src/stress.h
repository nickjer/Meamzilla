/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifndef MEAMZ_STRESS_H
#define MEAMZ_STRESS_H

#include "vect.h"
#include <istream>

namespace MEAMZ_NS {

class Stress
{
public:
  double xx, yy, zz, xy, yz, zx;

  Stress();
  Stress(const Stress&);
  Stress(double, double, double, double, double, double);
  virtual ~Stress();

  double& operator[](const int);                                // Return xx, yy, ...

  Stress& operator=(const double);                              // Set stress to value of double
  Stress& operator+=(const Stress&);                            // Compound assignment using addition with stress
  Stress& operator-=(const Stress&);                            // Compound assignment using subtraction with stress
  friend Stress operator&(const Vect&, const Vect&);            // Cross product

  friend std::istream& operator>>(std::istream&, Stress&);      // Input stream for stress

protected:

private:

};

/* ---------------------------------------------------------------------- */
inline
Stress::Stress() : xx(0), yy(0), zz(0), xy(0), yz(0), zx(0)
{
}

/* ---------------------------------------------------------------------- */
inline
Stress::Stress(const Stress& stress) : xx(stress.xx), yy(stress.yy), zz(stress.zz),
                                       xy(stress.xy), yz(stress.yz), zx(stress.zx)
{
}

/* ---------------------------------------------------------------------- */
inline
Stress::Stress(double _xx, double _yy, double _zz, double _xy, double _yz, double _zx)
  : xx(_xx), yy(_yy), zz(_zz), xy(_xy), yz(_yz), zx(_zx)
{
}

/* ---------------------------------------------------------------------- */
inline
Stress::~Stress()
{
}

/* ----------------------------------------------------------------------
   Return xx, yy, ... using brackets
------------------------------------------------------------------------- */
inline
double& Stress::operator[](const int idx)
{
  if (idx==0) return xx;
  else if (idx==1) return yy;
  else if (idx==2) return zz;
  else if (idx==3) return xy;
  else if (idx==4) return yz;
  else return zx;
}

/* ----------------------------------------------------------------------
   Set stress to value of double
------------------------------------------------------------------------- */
inline
Stress& Stress::operator=(const double rhs)
{
  xx = rhs;
  yy = rhs;
  zz = rhs;
  xy = rhs;
  yz = rhs;
  zx = rhs;
  return *this;
}

/* ----------------------------------------------------------------------
   Compound assignment using addition with stress
------------------------------------------------------------------------- */
inline
Stress& Stress::operator+=(const Stress &rhs)
{
  xx += rhs.xx;
  yy += rhs.yy;
  zz += rhs.zz;
  xy += rhs.xy;
  yz += rhs.yz;
  zx += rhs.zx;
  return *this;
}

/* ----------------------------------------------------------------------
   Compound assignment using subtraction with stress
------------------------------------------------------------------------- */
inline
Stress& Stress::operator-=(const Stress &rhs)
{
  xx -= rhs.xx;
  yy -= rhs.yy;
  zz -= rhs.zz;
  xy -= rhs.xy;
  yz -= rhs.yz;
  zx -= rhs.zx;
  return *this;
}

/* ----------------------------------------------------------------------
   Formation of stress tensor from vector product
------------------------------------------------------------------------- */
inline
Stress operator&(const Vect& lhs, const Vect& rhs)
{
  return Stress(lhs.x*rhs.x,lhs.y*rhs.y,lhs.z*rhs.z,lhs.x*rhs.y,lhs.y*rhs.z,lhs.z*rhs.x);
}

/* ----------------------------------------------------------------------
   Input stream for vect
------------------------------------------------------------------------- */
inline
std::istream& operator>>(std::istream& str, Stress& stress)
{
  return str >> stress.xx >> stress.yy >> stress.zz
             >> stress.xy >> stress.yz >> stress.zx;
}

}  // namespace MEAMZ_NS

#endif // MEAMZ_STRESS_H
