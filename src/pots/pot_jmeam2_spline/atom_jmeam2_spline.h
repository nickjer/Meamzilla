#ifdef ATOM_CLASS

AtomStyle(jmeam2/spline,AtomJMEAM2Spline)

#else

#ifndef MEAMZ_ATOM_JMEAM2_SPLINE_H
#define MEAMZ_ATOM_JMEAM2_SPLINE_H

#include "../pot_jmeam2/atom_jmeam2.h"

namespace MEAMZ_NS
{

class AtomJMEAM2Spline : public AtomJMEAM2
{
public:
  AtomJMEAM2Spline();
  virtual ~AtomJMEAM2Spline();

  virtual AtomJMEAM2Spline* clone() const;  // Copy polymorphic objects

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_JMEAM2_SPLINE_H
#endif // ATOM_CLASS
