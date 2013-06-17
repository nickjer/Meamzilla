#ifndef MEAMZ_ATOM_VEC_H
#define MEAMZ_ATOM_VEC_H

#include "pointers.h"
#include "atom.h"
#include <vector>

namespace MEAMZ_NS
{

class AtomVec : protected Pointers
{
public:
  int natoms;
  std::vector<Atom *> atoms;

  AtomVec(class Meamzilla *);
  virtual ~AtomVec();

  void init();              // Initialize atom vector

  bool operator<(const AtomVec& av) const {return nclusters_ < av.nclusters_;}

protected:

private:
  int npairs_;       // Total number of pairs this vector of atoms contains
  int ntriplets_;    // Total number of triplets this vector of atoms contains
  int nclusters_;    // Total number of pairs or triplets depending on what the potential sorts the atoms with

  void push_back(Atom *);   // Adds atom to list of atoms in this object
  void erase(int idx);      // Erases atom from list of atoms in this object

};

} // namespace MEAMZ_NS

#endif // MEAMZ_ATOM_VEC_H
