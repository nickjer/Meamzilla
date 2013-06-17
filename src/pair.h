#ifndef MEAMZ_PAIR_H
#define MEAMZ_PAIR_H

#include "vect.h"
#include "potential.h"

namespace MEAMZ_NS
{

class Pair
{
public:

  class Atom *neigh;    // Index of neigh atom corresponding to list of atoms in cell

  Vect dist;      // 3D radial "unit" vector between pair
  double r;       // radial distance
  double invr;    // inverse radial distance

  Pair();
  virtual ~Pair();

  virtual Pair* clone() const = 0;                  // Copy polymorphic objects

  void setup_pair(Atom *, Vect);                    // Setup pair with user specified conditions
  virtual bool check_pair(Atom *, Potential *) = 0; // Check that this pair is within boundaries of potentials

protected:

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_PAIR_H
