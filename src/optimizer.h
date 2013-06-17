#ifndef MEAMZ_OPTIMIZER_H
#define MEAMZ_OPTIMIZER_H

#include "pointers.h"
#include "potential.h"
#include "pot_list.h"
#include <ostream>
#include <vector>

namespace MEAMZ_NS
{

template<typename T>
using Vector = std::vector<T>;

class Optimizer : protected Pointers
{
public:
  Optimizer(class Meamzilla *, std::ostream * = nullptr);
  virtual ~Optimizer();

  virtual void init() = 0;              // Initialize the optimizer
  virtual void compute(PotList&) = 0;   // Use the optimizer on a potential list

protected:
  std::ostream *out_;                   // Output stream
  int max_steps_;                       // Maximum number of steps to perform in optimization

private:

};

} // namespace MEAMZ_NS

#endif // MEAMZ_OPTIMIZER_H
