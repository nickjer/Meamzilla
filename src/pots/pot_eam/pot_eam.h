#ifdef POTENTIAL_CLASS

PotentialStyle(eam,PotEAM)

#else

#ifndef MEAMZ_POT_EAM_H
#define MEAMZ_POT_EAM_H

#include "../pot_pair/pot_pair.h"

namespace MEAMZ_NS {

class PotEAM : public PotPair
{
public:
  PotEAM(class Meamzilla *, int);
  virtual ~PotEAM();

  virtual PotEAM* clone() const;                      // Copy polymorphic objects

  virtual void init();                                // Initialize potential

  virtual void check_if_similar(Potential *) const;   // Check that a potential is similar enough to this
                                                      // potential for pairs/triplets to have identical info

  virtual void compute_trap(const Comm&, int);                    // Trap non-root procs in infinite loop
  virtual double compute(const Comm&, ErrorVec * = nullptr);      // Compute error sum/vector
  virtual int rescale(const Comm&, std::ostream *, int = 0);      // Rescale potential

protected:
  int embed_extrap_;   // Flag of whether or not we linearly extrapolate out embedding fn in force computation

  virtual void compute_densities(const Comm&, Vector<double> *);  // Compute densities for each atom
  virtual void write_punish(const ErrorVec&, String, int) const;  // Write punishment/constraint errors (none for this pot type)
  virtual void write_extras(const Comm&, String, int);            // Write out extra data (densities for this pot)
  virtual void write_densities(const Comm&, String, int);         // Write out densities

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_POT_EAM_H
#endif // POTENTIAL_CLASS
