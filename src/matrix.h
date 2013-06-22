/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   A fully contiguous Matrix defined as matrix(row,col), where we
   keep the columns contiguous in memory for Lapack usage
------------------------------------------------------------------------- */

#ifndef MEAMZ_MATRIX_H
#define MEAMZ_MATRIX_H

#include <cmath>
#include <algorithm>
#include <cstring>

namespace MEAMZ_NS {

template <class T>
class Matrix {
public:
  int m_stride;
  int m_size;
  T *m_data;

  Matrix();
  Matrix(const Matrix&);
  Matrix(int, int);
  Matrix(int, int, int );
  ~Matrix();

  void resize(int, int);
  T& operator()(int, int);
  T operator()(int, int) const;
  T& operator=(const T&);

protected:

private:

};

/* ---------------------------------------------------------------------- */
template <class T>
inline
Matrix<T>::Matrix() : m_stride(0), m_size(0), m_data(nullptr)
{
}

/* ---------------------------------------------------------------------- */
template <class T>
inline
Matrix<T>::Matrix(const Matrix& other) : m_stride(other.m_stride), m_size(other.m_size), m_data(nullptr)
{
  m_data = new T [other.m_size];
  std::copy(other.m_data, other.m_data + other.m_size, m_data);
}

/* ---------------------------------------------------------------------- */
template <class T>
inline
Matrix<T>::Matrix(int dimY, int dimX) : m_stride(dimY), m_size(dimX * dimY), m_data(new T[dimX * dimY])
{
}

/* ---------------------------------------------------------------------- */
template <class T>
inline
Matrix<T>::Matrix(int dimY, int dimX, int value) : m_stride(dimY), m_size(dimX * dimY), m_data(new T[dimX * dimY])
{
  memset(m_data, value, sizeof(T)*m_size);
}

/* ---------------------------------------------------------------------- */
template <class T>
inline
Matrix<T>::~Matrix()
{
  delete [] m_data;
}

/* ----------------------------------------------------------------------
   Resize the matrix (needs optimization work, but this isn't
   called all that often)
------------------------------------------------------------------------- */
template <class T>
inline
void Matrix<T>::resize(int dimY, int dimX)
{
  m_stride = dimY;
  m_size = dimX * dimY;
  if (m_data) delete [] m_data;
  m_data = new T [m_size];
}

/* ----------------------------------------------------------------------
   Set vector to value of double
------------------------------------------------------------------------- */
template <class T>
inline
T& Matrix<T>::operator()(int row, int col)
{
  return *(m_data + row + m_stride * col);
}

/* ----------------------------------------------------------------------
   Read value
------------------------------------------------------------------------- */
template <class T>
inline
T Matrix<T>::operator()(int row, int col) const
{
  return *(m_data + row + m_stride * col);
}

/* ----------------------------------------------------------------------
   Assignment operator
------------------------------------------------------------------------- */
template <class T>
inline
T& Matrix<T>::operator=(const T& other)
{
  // protect against invalid self-assignment
  if (this != &other) {
      // 1: allocate new memory and copy the elements
      T *new_data = new T [other.m_size];
      std::copy(other.m_data, other.m_data + other.m_size, new_data);

      // 2: deallocate old memory
      delete [] m_data;

      // 3: assign the new memory to the object
      m_data = new_data;
      m_size = other.m_size;
      m_stride = other.m_stride;
  }
  // by convention, always return *this
  return *this;
}

} // MEAMZ_NS

#endif // MEAMZ_MATRIX_H
