#include "opt_genalg.h"
#include "opt_powell.h"
#include "input.h"
#include "random.h"
#include "timer.h"

#include <iomanip>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

OptGenAlg::OptGenAlg(class Meamzilla *mmz, std::ostream *out) : Optimizer(mmz, out), pop_size_(10),
                                                                cross_rate_(0.1), mut_rate_(0.0),
                                                                init_scale_(10.0), fit_rate_(0.2), rescale_rate_(0.1),
                                                                order_breed_(0), gen_save_(0), num_powell_(0)
{
  max_steps_ = 1000;
}

/* ---------------------------------------------------------------------- */

OptGenAlg::~OptGenAlg()
{
  //dtor
}

/* ----------------------------------------------------------------------
   initialize optimizer
------------------------------------------------------------------------- */

void OptGenAlg::init()
{
  if (mmz->input) {
    mmz->input->parse("pop_size", 0, pop_size_);
    mmz->input->parse("cross_rate", 0, cross_rate_);
    mmz->input->parse("mut_rate", 0, mut_rate_);
    mmz->input->parse("init_scale", 0, init_scale_);
    mmz->input->parse("num_powell", 0, num_powell_);
    mmz->input->parse("fit_rate", 0, fit_rate_);
    mmz->input->parse("rescale_rate", 0, rescale_rate_);
    mmz->input->parse("order_breed", 0, order_breed_);
    mmz->input->parse("gen_save", 0, gen_save_);
  }

  if (pop_size_ < 3) {
    if (mmz->universe->is_root())
      std::cerr << "ERROR: Need a population size larger than 2 for genalg" << std::endl;
    mmz->universe->abort_all();
  }

  return;
}

/* ----------------------------------------------------------------------
   use the optimizer on a potential list
------------------------------------------------------------------------- */

void OptGenAlg::compute(PotList& pot_list)
{
  if (mmz->universe->is_root() && out_) {
    *out_ << "Beginning genetic algorithm minimization..." << std::endl << std::endl;
  }
  mmz->universe->barrier();

  // Treat entire potlist as population and minimize it
  mmz->random->reset_seed();    // reset the seed to original seed for reproducibility
  pop_size_ = std::max(pop_size_, pot_list.get_npots());  // set population size to larger number

  Potential *tmp_pot = pot_list.pot_template()->clone();
  tmp_pot->compute_trap(mmz->universe->comm_inter, 1);    // trap procs that aren't root of comm_inter
  if (mmz->universe->comm_inter.is_root()) {
    pot_list = compute_potlist(pot_list, max_steps_);
    tmp_pot->compute_trap(mmz->universe->comm_inter, 0);    // untrap procs
  }
  delete tmp_pot;

  // Only save the top number of potentials the user specified and dump the rest
  if (gen_save_ > 0) pot_list.resize(gen_save_);

  pot_list.write_pots(pot_list.get_end_file() + ".genalg"); // this is also a barrier for procs

  if (mmz->universe->is_root() && out_)
    *out_ << std::endl << "Finished genetic algorithm minimization..." << std::endl << std::endl;
  mmz->universe->barrier();

  return;
}

/* ----------------------------------------------------------------------
   optimization loop
------------------------------------------------------------------------- */

PotList OptGenAlg::compute_potlist(const PotList& pot_list, int max_steps) const
{
  // Some local variables to make my life easier
  PotList opt_pot_list = pot_list;                    // Optimized population list
  Comm& comm_all = mmz->universe->comm_all;           // Communicator for all procs
  Comm& comm_inter = mmz->universe->comm_inter;       // Communicator for procs in a group
  Comm& comm_intra = mmz->universe->comm_intra;       // Communicator for root procs in each group
  int group = mmz->universe->get_group();             // Group number
  int ncoeff = pot_list.pot_template()->get_ncoeff(); // Number of coefficients to be optimized (taken from first pot in list)

  int loopagain;
  Vector<double> F1var(max_steps+NEPS,1.e300);   // Backlog of old F1 values
  Vector<double> F2var(max_steps+NEPS,1.e300);   // Backlog of old F2 values

  // Add extra potentials to pot_list (modified first potential)
  if (pot_list.get_npots() < pop_size_) {
    int remainder = pop_size_ - pot_list.get_npots();

    // Modify first potential and add it to list until list is pop_size_ length
    for (int i=0; i<remainder; ++i) {
      Potential *pot = pot_list.pot_template()->clone();

      // Modify each coefficient by gaussian distribution with scale = coeff_value * init_scale
      for (int c=0; c<ncoeff; ++c) {
        double height = mmz->random->norm_dist() * pot->coeff(c) * init_scale_;
        double width = std::abs(mmz->random->norm_dist());
        pot->modify_coeff(c, height, width);
      }

      opt_pot_list.push_back(pot);
    }
  }

  // Compute error sums on each potential in list and output it
  if (comm_all.is_root() && out_)
    *out_ << "Initializing population with population size = " << pop_size_ << std::endl;

  int pop_size = opt_pot_list.get_group_npots();  // Local population size seen by a single group
  for (int i=0; i<pop_size; ++i) {
    double error_sum = opt_pot_list[i]->compute(comm_inter);
    if (out_) {
      int pot_idx = opt_pot_list[i]->get_global_idx();
      *out_ << "G#" << group << " P#" << pot_idx << " - ";
      *out_ << std::scientific << std::setprecision(8) << error_sum << std::endl;
    }
  }

  opt_pot_list.sort();  // potentials are communicated to all procs along with error_sum, then sorted

  // Begin looping through genetic algorithm steps
  int count = 0;
  do {
    Timer timer_genalg;  // setup timer and start it
    auto ncalls_begin = Potential::get_ncalls();  // Get number of compute calls for potential

    if (comm_all.is_root() && out_) {
      *out_ << std::endl;
      *out_ << "Step: " << count+1 << std::endl;
      *out_ << "-------------------------------" << std::endl;
      *out_ << "Breeding ( 0+0 1+1 ";
    }

    comm_intra.barrier();

    PotList tmp_pot_list = opt_pot_list;  // temporary potlist for children to go

    // Make list of potential indices that we use later on when breeding
    int parent_idx = 0;
    Vector<int> pot_idx_list(pop_size_-2);
    for (int i=0; i<pop_size_-2; ++i) pot_idx_list[i] = i+1;

    for (int i=0; i<pop_size_; ++i) {
      Potential*& new_pot = tmp_pot_list.global_at(i);

      // Initialize this potential to current potential
      delete new_pot;
      new_pot = opt_pot_list.global_at(i)->clone();
      new_pot->set_global_idx(i);

      // Copy over the best 2 parents
      if (i==0 || i==1) continue;

      // Parent to breed with a potential
      if (i>2 && mmz->random->unif_dist() < fit_rate_) ++parent_idx;

      // Potential to be bred with parent
      int parent2_idx;
      if (order_breed_) {  // breed potentials in order of quality with parent_idx
        parent2_idx = i-1;
      } else {  // randomly breed potential with parent_idx
        int temp_idx;
        do {
          temp_idx = mmz->random->unif_dist() * pot_idx_list.size();
          parent2_idx = pot_idx_list[temp_idx];
          if (pot_idx_list.size() == 1 && parent2_idx == parent_idx) ++parent_idx;
        } while (parent2_idx == parent_idx);
        pot_idx_list.erase(pot_idx_list.begin() + temp_idx); // We already used this pot, so erase it
      }

      // Crossover operation for two parents
      int parent_switch = mmz->random->unif_dist() * 2; // 0 or 1
      int success = 0;
      for (int c=0; c<ncoeff; ++c) {
        // Decide if switch parent to seed the DNA of child
        if (mmz->random->unif_dist() < cross_rate_) {
          parent_switch = (parent_switch + 1) % 2;
          success = 1;
        }
        if (c == ncoeff-1 && success == 0) parent_switch = (parent_switch + 1) % 2;   // Force at least one crossover
        int pot_idx = (parent_switch==0) ? parent_idx : parent2_idx;  // either use parent or parent2
        new_pot->coeff(c) = opt_pot_list.global_at(pot_idx)->coeff(c);

        // Mutation (take coeff from random potential in population)
        if (mmz->random->unif_dist() < mut_rate_) {
          int rand_pot_idx;
          do {
            rand_pot_idx = mmz->random->unif_dist() * pop_size_;
          } while (rand_pot_idx == parent_idx || rand_pot_idx == parent2_idx); // ignore parent potentials
          new_pot->coeff(c) = opt_pot_list.global_at(rand_pot_idx)->coeff(c);
          //new_pot->coeff(c) += new_pot->coeff(c) * 0.1 * mmz->random->norm_dist();
        }
      }

      if (comm_all.is_root() && out_)
        *out_ << parent_idx << "+" << parent2_idx << " ";
    }

    if (comm_all.is_root() && out_)
      *out_ << ")" << std::endl << std::endl;

    opt_pot_list = tmp_pot_list;

    for (int i=0; i<pop_size; ++i) {
      int did_rescale = 0;
      if (mmz->random->unif_dist() < rescale_rate_ && i>1) {
        opt_pot_list[i]->rescale(comm_inter, out_, 1);
        did_rescale = 1;
      }
      opt_pot_list[i]->compute(comm_inter);

      // Create local optimizer
      OptPowell optimizer_pow(mmz);
      optimizer_pow.init();
      Potential *pot = optimizer_pow.compute_pot(opt_pot_list[i], num_powell_);
      if (pot->get_error() <= opt_pot_list[i]->get_error()) {
        delete opt_pot_list[i];
        opt_pot_list[i] = pot->clone();
      }
      delete pot;

      double error_sum = opt_pot_list[i]->get_error();

      if (out_) {
        int pot_idx = opt_pot_list[i]->get_global_idx();
        *out_ << "G#" << group << " P#" << pot_idx << " - ";
        *out_ << std::scientific << std::setprecision(8) << error_sum;
        if (did_rescale)
          *out_ << "  (rescaled)";
        *out_ << std::endl;
      }
    }

    timer_genalg.stop();  // stop the timer
    if (out_) {
      auto ncalls_end = Potential::get_ncalls();
      double time = timer_genalg.get_time();  // get the recorded time from start to stop
      *out_ << "G#" << group << " P#**" << " - ";
      *out_ << "TTT - calls = " << ncalls_end - ncalls_begin
            << ", total time = " << std::fixed << time
            << " sec, calls/sec = " << (ncalls_end - ncalls_begin)*1.0/time << std::endl;
    }

    opt_pot_list.sort();  // potentials are communicated to all procs along with error_sum, then sorted

    opt_pot_list.write_pots(mmz->potlist->get_tmp_file() + ".genalg");

    comm_intra.barrier();

    // Drop out of loop if last few top two F's aren't changing much
    loopagain = 0;
    F1var[count + NEPS] = opt_pot_list.global_at(0)->get_error(); // record best F in backlog of F's
    F2var[count + NEPS] = opt_pot_list.global_at(1)->get_error(); // record 2nd best F in backlog of F's
    double F1 = opt_pot_list.global_at(0)->get_error();
    for (int i = 1; i <= NEPS; ++i)
      if (std::abs(F1 - F1var[count + NEPS - i]) > (EPS * F1)) loopagain = 1;
    double F2 = opt_pot_list.global_at(1)->get_error();
    for (int i = 1; i <= NEPS; ++i)
      if (std::abs(F2 - F2var[count + NEPS - i]) > (EPS * F2)) loopagain = 1;

    ++count;
  } while (count<max_steps && loopagain);

  if (comm_all.is_root() && out_) {
    *out_ << std::endl;
    *out_ << "Final Potentials: " << std::endl;
    *out_ << "-------------------------------" << std::endl;
    for (int i=0; i<pop_size_; ++i) {
      double error_sum = opt_pot_list.global_at(i)->get_error();
      *out_ << "P#" << i << " - ";
      *out_ << std::scientific << std::setprecision(8) << error_sum << std::endl;
    }
  }

  return opt_pot_list;
}
