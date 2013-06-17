#ifndef MEAMZ_OPT_LIST_H
#define MEAMZ_OPT_LIST_H

#include "pointers.h"
#include "optimizer.h"

namespace MEAMZ_NS
{

class OptList : protected Pointers
{
public:
  OptList(class Meamzilla *);
  virtual ~OptList();

  void init();                  // Initialize this class
  void compute(PotList&);       // Compute list of optimizers using list of potentials

protected:

private:
  int nopts_;                       // Number of potentials
  std::vector<Optimizer *> opts_;   // List of optimizers

};

} // namespace MEAMZ_NS

#endif // MEAMZ_OPT_LIST_H
