#ifdef ATOM_CLASS

AtomStyle(pair,AtomPair)

#else

#ifndef MEAMZ_ATOM_PAIR_H
#define MEAMZ_ATOM_PAIR_H

#include "../../atom.h"

namespace MEAMZ_NS
{

class AtomPair : public Atom
{
public:
  AtomPair();
  virtual ~AtomPair();

  virtual AtomPair* clone() const;  // Copy polymorphic objects

  virtual int cmp_value() const;    // Value used when sorting atoms

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_PAIR_H
#endif // ATOM_CLASS
