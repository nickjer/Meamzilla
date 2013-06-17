#ifdef ATOM_CLASS

AtomStyle(meam/spline,AtomMEAMSpline)

#else

#ifndef MEAMZ_ATOM_MEAM_SPLINE_H
#define MEAMZ_ATOM_MEAM_SPLINE_H

#include "../pot_meam/atom_meam.h"

namespace MEAMZ_NS
{

class AtomMEAMSpline : public AtomMEAM
{
public:
  AtomMEAMSpline();
  virtual ~AtomMEAMSpline();

  virtual AtomMEAMSpline* clone() const;  // Copy polymorphic objects

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_MEAM_SPLINE_H
#endif // ATOM_CLASS
