#ifdef ATOM_CLASS

AtomStyle(eam,AtomEAM)

#else

#ifndef MEAMZ_ATOM_EAM_H
#define MEAMZ_ATOM_EAM_H

#include "../pot_pair/atom_pair.h"

namespace MEAMZ_NS
{

class AtomEAM : public AtomPair
{
public:
  int F_idx;              // index of embedding function for this atom (equal to atom type)

  AtomEAM();
  virtual ~AtomEAM();

  virtual AtomEAM* clone() const;  // Copy polymorphic objects

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_EAM_H
#endif // ATOM_CLASS
