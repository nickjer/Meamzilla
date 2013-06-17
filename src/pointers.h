#ifndef MEAMZ_POINTERS_H
#define MEAMZ_POINTERS_H

#include "meamzilla.h"

namespace MEAMZ_NS
{

class Pointers
{
public:
  Pointers(Meamzilla *ptr) : mmz(ptr) {}
  virtual ~Pointers() {}

protected:
  Meamzilla *mmz;

private:
};

} // namespace MEAMZ_NS

#endif // MEAMZ_POINTERS_H
