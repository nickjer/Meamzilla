#ifdef ATOM_CLASS

AtomStyle(jmeam,AtomJMEAM)

#else

#ifndef MEAMZ_ATOM_JMEAM_H
#define MEAMZ_ATOM_JMEAM_H

#include "../pot_meam/atom_meam.h"

namespace MEAMZ_NS
{

class AtomJMEAM : public AtomMEAM
{
public:
  AtomJMEAM();
  virtual ~AtomJMEAM();

  virtual AtomJMEAM* clone() const;  // Copy polymorphic objects

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_JMEAM_H
#endif // ATOM_CLASS
