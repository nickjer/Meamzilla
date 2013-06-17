#ifndef MEAMZ_ERROR_H
#define MEAMZ_ERROR_H

#include <string>
#include <stdexcept>

namespace MEAMZ_NS
{

class Error : public std::runtime_error
{
public:
  Error(const std::string& msg) : std::runtime_error(msg) {}
};

} // namespace MEAMZ_NS

#endif // MEAMZ_ERROR_H
