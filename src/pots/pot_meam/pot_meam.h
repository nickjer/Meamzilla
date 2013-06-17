#ifdef POTENTIAL_CLASS

PotentialStyle(meam,PotMEAM)

#else

#ifndef MEAMZ_POT_MEAM_H
#define MEAMZ_POT_MEAM_H

#include "../pot_eam/pot_eam.h"

namespace MEAMZ_NS {

class PotMEAM : public PotEAM
{
public:
  PotMEAM(class Meamzilla *, int);
  virtual ~PotMEAM();

  virtual PotMEAM* clone() const;                     // Copy polymorphic objects

  virtual void check_if_similar(Potential *) const;   // Check that a potential is similar enough to this
                                                      // potential for pairs/triplets to have identical info

  virtual double compute(const Comm&, ErrorVec * = nullptr);      // Compute error sum/vector
  virtual int rescale(const Comm&, std::ostream *, int = 0);      // Rescale potential

protected:
  virtual void compute_densities(const Comm&, Vector<double> *);  // Compute densities for each atom
  virtual void write_punish(const ErrorVec&, String, int) const;  // Write punishment/constraint errors (none for this pot type)

  virtual int rescale_3body(const Comm&, std::ostream *, int);    // Rescale MEAM part of the potential: f*f*g

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_POT_MEAM_H
#endif // POTENTIAL_CLASS
