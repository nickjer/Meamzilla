#ifndef MEAMZ_RANDOM_H
#define MEAMZ_RANDOM_H

#include "pointers.h"
#include <functional>
#include <random>

namespace MEAMZ_NS
{

class Random : protected Pointers
{
public:
  Random(class Meamzilla *);
  virtual ~Random();

  std::function<double()> unif_dist;    // Return a uniform distribution [0:1)
  std::function<double()> norm_dist;    // Return a normal distribution with mean=0 and std=1

  void init();
  void reset_seed();                    // Reset the seed to user specified seed

protected:

private:
  unsigned long seed_;  // seed

  std::mt19937 rng_;    // Engine used (Mersenne Twister)

  // Distributions used
  std::uniform_real_distribution<double> ureal_dist_;
  std::normal_distribution<double> normal_dist_;

};

} // namespace MEAMZ_NS

#endif // MEAMZ_RANDOM_H
