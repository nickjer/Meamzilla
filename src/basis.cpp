/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "basis.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Basis::Basis() : ncoeff_(0), basis_type_("")
{
  //ctor
}

/* ---------------------------------------------------------------------- */

Basis::~Basis()
{
  //dtor
}

/* ----------------------------------------------------------------------
   get basis type for this basis fn
------------------------------------------------------------------------- */

String Basis::get_basis_type() const
{
  return basis_type_;
}

/* ----------------------------------------------------------------------
   get number of adjustable coefficients for this basis fn
------------------------------------------------------------------------- */

int Basis::get_ncoeff() const
{
  return ncoeff_;
}
