#include "random.h"
#include "input.h"

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Random::Random(class Meamzilla *mmz) : Pointers(mmz), seed_(0), ureal_dist_(0,1), normal_dist_(0,1)
{
  unif_dist = std::bind( ureal_dist_, std::ref(rng_) );
  norm_dist = std::bind( normal_dist_, std::ref(rng_) );
}

/* ---------------------------------------------------------------------- */

Random::~Random()
{
  //dtor
}

/* ----------------------------------------------------------------------
   initialize random number generator with user supplied seed
------------------------------------------------------------------------- */

void Random::init()
{
  if (mmz->input)
    mmz->input->parse("seed", 0, seed_);

  // NOTE: Each processor will have the same seed, if you'd like
  // different seeds for different groups of processors then I
  // suggest adding the group number to the seed before seeding
  // the generator. But this will make reproducing results
  // impossible unless user uses same number of groups second
  // time around.

  rng_.seed(seed_);   // seed the generator
}

/* ----------------------------------------------------------------------
   reset seed to user specified seed
------------------------------------------------------------------------- */

void Random::reset_seed()
{
  rng_.seed(seed_);   // seed the generator
}
