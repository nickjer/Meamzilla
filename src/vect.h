#ifndef MEAMZ_VECT_H
#define MEAMZ_VECT_H

#include <cmath>
#include <istream>

namespace MEAMZ_NS {

class Vect {
public:
  double x, y, z;

  Vect();
  Vect(const Vect&);
  Vect(double, double, double);
  ~Vect();

  double& operator[](const int);                                // Return x, y, or z

  Vect& operator=(const double);                                // Set vector to value of double
  Vect& operator+=(const Vect&);                                // Compound assignment using addition with vector
  Vect& operator-=(const Vect&);                                // Compound assignment using subtraction with vector
  Vect& operator*=(const double);                               // Compound assignment using multiplication with double
  friend double operator*(const Vect&, const Vect&);            // Dot product
  friend Vect operator+(const Vect&, const Vect&);              // Add vectors
  friend Vect operator-(const Vect&, const Vect&);              // Subtract vectors
  friend Vect operator^(const Vect&, const Vect&);              // Cross product
  friend Vect operator*(const Vect&, const double);             // Multiply vector by double
  friend Vect operator*(const double, const Vect&);             // Multiply vector by double
  friend Vect operator/(const Vect&, const double);             // Divide vector by double
  friend double mag(const Vect&);                               // Magnitude of vecctor
  friend std::istream& operator>>(std::istream&, Vect&);        // Input stream for vect

protected:

private:

};

/* ---------------------------------------------------------------------- */
inline
Vect::Vect() : x(0), y(0), z(0)
{
}

/* ---------------------------------------------------------------------- */
inline
Vect::Vect(const Vect& vect) : x(vect.x), y(vect.y), z(vect.z)
{
}

/* ---------------------------------------------------------------------- */
inline
Vect::Vect(double _x, double _y, double _z) : x(_x), y(_y), z(_z)
{
}

/* ---------------------------------------------------------------------- */
inline
Vect::~Vect()
{
}

/* ----------------------------------------------------------------------
   Return x, y, or z using brackets
------------------------------------------------------------------------- */
inline
double& Vect::operator[](const int idx)
{
  if (idx==0) return x;
  else if (idx==1) return y;
  else return z;
}

/* ----------------------------------------------------------------------
   Set vector to value of double
------------------------------------------------------------------------- */
inline
Vect& Vect::operator=(const double rhs)
{
  x = rhs;
  y = rhs;
  z = rhs;
  return *this;
}

/* ----------------------------------------------------------------------
   Compound assignment using addition with vector
------------------------------------------------------------------------- */
inline
Vect& Vect::operator+=(const Vect &rhs)
{
  x += rhs.x;
  y += rhs.y;
  z += rhs.z;
  return *this;
}

/* ----------------------------------------------------------------------
   Compound assignment using subtraction with vector
------------------------------------------------------------------------- */
inline
Vect& Vect::operator-=(const Vect &rhs)
{
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;
  return *this;
}

/* ----------------------------------------------------------------------
   Compound assignment using multiplication with double
------------------------------------------------------------------------- */
inline
Vect& Vect::operator*=(const double rhs)
{
  x *= rhs;
  y *= rhs;
  z *= rhs;
  return *this;
}

/* ----------------------------------------------------------------------
   Dot Product
------------------------------------------------------------------------- */
inline
double operator*(const Vect& lhs, const Vect& rhs)
{
  return lhs.x*rhs.x+lhs.y*rhs.y+lhs.z*rhs.z;
}

/* ----------------------------------------------------------------------
   Add vectors
------------------------------------------------------------------------- */
inline
Vect operator+(const Vect& lhs, const Vect& rhs)
{
  return Vect(lhs.x+rhs.x,lhs.y+rhs.y,lhs.z+rhs.z);
}

/* ----------------------------------------------------------------------
   Subtract vectors
------------------------------------------------------------------------- */
inline
Vect operator-(const Vect& lhs, const Vect& rhs)
{
  return Vect(lhs.x-rhs.x,lhs.y-rhs.y,lhs.z-rhs.z);
}

/* ----------------------------------------------------------------------
   Cross Product
------------------------------------------------------------------------- */
inline
Vect operator^(const Vect& lhs, const Vect& rhs)
{
  return Vect(lhs.y*rhs.z-lhs.z*rhs.y,lhs.z*rhs.x-lhs.x*rhs.z,lhs.x*rhs.y-lhs.y*rhs.x);
}

/* ----------------------------------------------------------------------
   Multiply vector by double
------------------------------------------------------------------------- */
inline
Vect operator*(const Vect& lhs, const double rhs)
{
  return Vect(lhs.x*rhs,lhs.y*rhs,lhs.z*rhs);
}

/* ----------------------------------------------------------------------
   Multiply vector by double
------------------------------------------------------------------------- */
inline
Vect operator*(const double lhs, const Vect& rhs)
{
  return rhs*lhs;
}

/* ----------------------------------------------------------------------
   Divide vector by double
------------------------------------------------------------------------- */
inline
Vect operator/(const Vect& lhs, const double rhs)
{
  return Vect(lhs.x/rhs,lhs.y/rhs,lhs.z/rhs);
}

/* ----------------------------------------------------------------------
   Magnitude of vector
------------------------------------------------------------------------- */
inline
double mag(const Vect& v)
{
  return std::sqrt(v*v);
}

/* ----------------------------------------------------------------------
   Input stream for vect
------------------------------------------------------------------------- */
inline
std::istream& operator>>(std::istream& str, Vect& vect)
{
  return str >> vect.x >> vect.y >> vect.z;
}

}  // namespace MEAMZ_NS

#endif // MEAMZ_VECT_H
