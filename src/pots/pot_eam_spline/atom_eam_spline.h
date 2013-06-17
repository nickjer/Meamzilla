#ifdef ATOM_CLASS

AtomStyle(eam/spline,AtomEAMSpline)

#else

#ifndef MEAMZ_ATOM_EAM_SPLINE_H
#define MEAMZ_ATOM_EAM_SPLINE_H

#include "../pot_eam/atom_eam.h"

namespace MEAMZ_NS
{

class AtomEAMSpline : public AtomEAM
{
public:
  AtomEAMSpline();
  virtual ~AtomEAMSpline();

  virtual AtomEAMSpline* clone() const;  // Copy polymorphic objects

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_EAM_SPLINE_H
#endif // ATOM_CLASS
