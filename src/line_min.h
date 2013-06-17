#ifndef MEAMZ_LINE_MIN_H
#define MEAMZ_LINE_MIN_H

#include <functional>

namespace MEAMZ_NS
{

class LineMin
{
public:
  LineMin();
  virtual ~LineMin();

  double get_xmin1() const;
  double get_xmin2() const;
  double get_Fmin1() const;
  double get_Fmin2() const;

  void set_func(std::function<double(double)>&);

  virtual double find_minima() = 0;   // Find the minima along line, return error_sum

protected:
  std::function<double(double)> eval_;

  double xa_, xb_, xc_;     // Unitless values describing the proportional distance along the direction vector
  double Fa_, Fb_, Fc_;     // Corresponding function values at these points along line

  double xmin1_, xmin2_;    // Unitless values describing the distance along the direction vector for minima and 2nd minima
  double Fxmin1_, Fxmin2_;  // Corresponding fn values

  // Constants used in bracketing function
  static constexpr double CGOLD = 0.381966012;      // Golden ratio
  static constexpr unsigned int MAX_BRACKET = 100;  // Max # of bracketing steps

private:
  void get_bracket();       // Bracket the minimum (xa,xb,xc) where f(xa)>f(xb)<f(xc)

};

} // namespace MEAMZ_NS

#endif // MEAMZ_LINE_MIN_H
