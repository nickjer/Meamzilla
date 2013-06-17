#ifdef ATOM_CLASS

AtomStyle(pair/spline,AtomPairSpline)

#else

#ifndef MEAMZ_ATOM_PAIR_SPLINE_H
#define MEAMZ_ATOM_PAIR_SPLINE_H

#include "../pot_pair/atom_pair.h"

namespace MEAMZ_NS
{

class AtomPairSpline : public AtomPair
{
public:
  AtomPairSpline();
  virtual ~AtomPairSpline();

  virtual AtomPairSpline* clone() const;  // Copy polymorphic objects

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_PAIR_SPLINE_H
#endif // ATOM_CLASS
