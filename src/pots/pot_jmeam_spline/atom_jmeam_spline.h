#ifdef ATOM_CLASS

AtomStyle(jmeam/spline,AtomJMEAMSpline)

#else

#ifndef MEAMZ_ATOM_JMEAM_SPLINE_H
#define MEAMZ_ATOM_JMEAM_SPLINE_H

#include "../pot_jmeam/atom_jmeam.h"

namespace MEAMZ_NS
{

class AtomJMEAMSpline : public AtomJMEAM
{
public:
  AtomJMEAMSpline();
  virtual ~AtomJMEAMSpline();

  virtual AtomJMEAMSpline* clone() const;  // Copy polymorphic objects

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_JMEAM_SPLINE_H
#endif // ATOM_CLASS
