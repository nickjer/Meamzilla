/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "opt_simann.h"
#include "input.h"
#include "timer.h"
#include "random.h"

#include <cmath>
#include <iomanip>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

OptSimAnn::OptSimAnn(class Meamzilla *mmz, std::ostream *out) : Optimizer(mmz, out), anneal_temp_(0.0)
{
  max_steps_ = 1000;
}

/* ---------------------------------------------------------------------- */

OptSimAnn::~OptSimAnn()
{
  //dtor
}

/* ----------------------------------------------------------------------
   initialize optimizer
------------------------------------------------------------------------- */

void OptSimAnn::init()
{
  if (mmz->input) {
    mmz->input->parse("anneal_temp", 1, anneal_temp_);
  }

  return;
}

/* ----------------------------------------------------------------------
   use the optimizer on a potential list
------------------------------------------------------------------------- */

void OptSimAnn::compute(PotList& pot_list)
{
  if (anneal_temp_ <= 0.0) return;    // Check if anneal temp is set

  if (mmz->universe->is_root() && out_) {
    *out_ << "Beginning simulated annealing..." << std::endl << std::endl;
    print_simann(-1, 0, 0, 0, 0, 0, 0);
  }
  mmz->universe->barrier();

  // Loop through all the potentials this group sees
  int npots = pot_list.get_group_npots();
  for (int p=0; p<npots; ++p) {
    mmz->random->reset_seed();    // reset the seed to original seed for reproducibility
    Potential::reset_ncalls();    // reset the number of force computation calls for each potential
    pot_list[p]->compute_trap(mmz->universe->comm_inter, 1);  // trap procs that aren't root of comm_inter
    if (!mmz->universe->comm_inter.is_root()) break;
    Potential *pot = compute_pot(pot_list[p], max_steps_);
    pot_list[p]->compute_trap(mmz->universe->comm_inter, 0);  // untrap procs
    delete pot_list[p];
    pot_list[p] = pot;
  } // End loop over potentials

  pot_list.write_pots(pot_list.get_end_file() + ".simann"); // this is also a barrier for procs

  if (mmz->universe->is_root() && out_)
    *out_ << std::endl << "Finished simulated annealing..." << std::endl << std::endl;
  mmz->universe->barrier();

  return;
}

/* ----------------------------------------------------------------------
   optimization loop
------------------------------------------------------------------------- */

Potential* OptSimAnn::compute_pot(const Potential *main_pot, int max_steps) const
{
  // Some local variables to make my life easier
  Comm& comm_inter = mmz->universe->comm_inter;   // Communicator for procs in a group
  int group = mmz->universe->get_group();         // Group number
  int pot_idx = main_pot->get_global_idx();       // Global potential index in list of potentials for this pot
  int ncoeff = main_pot->get_ncoeff();            // Number of coefficients to be optimized

  Potential *pot = main_pot->clone();             // Clone main potential to local potential
  double F = pot->compute(comm_inter);            // Get force error

  if (max_steps == 0) return pot;

  Potential *opt_pot = pot->clone();              // Set up optimal potential
  double Fopt = F;                                // Set force error

  std::vector<double> v(ncoeff, 1.0);             // Scaling factor for each coefficient
  std::vector<int> naccept(ncoeff, 0);            // Number of times we accepted change in coefficient
  std::vector<double> Fvar(max_steps+5+NEPS,F);   // Backlog of old F values

  int k = 0;                  // Temperature step
  double T = anneal_temp_;    // Temperature
  print_simann(k, T, 0, F, Fopt, group, pot_idx);

  // Enter main simulated annealing loop
  int loopagain = 1;
  int rescale_me = 1; // turn on rescaling & allow it to shut off if reset local min to global min
  do {
    Timer timer_simann;  // setup timer and start it
    auto ncalls_begin = Potential::get_ncalls();  // Get number of compute calls for potential

    for (int m=0; m<3*ncoeff; ++m) {      // Number of steps at fixed temp
      for (int j=0; j<NSTEP; ++j) {       // Run through NSTEP optimizations of each coefficient
        for (int c=0; c<ncoeff; ++c) {

          Potential *trial_pot = pot->clone();   // setup trial_pot at each step

          double height = mmz->random->norm_dist() * v[c];
          double width = std::abs(mmz->random->norm_dist());
          trial_pot->modify_coeff(c, height, width);
          double F2 = trial_pot->compute(comm_inter); // compute error

          naccept[c] += metropolis(pot, trial_pot, &F, F2, T, opt_pot, &Fopt);

          delete trial_pot; // clean up
        }
      } // END OPTIMIZATION OF COEFFICIENT INNER AND OUTER LOOPS

      // Adjust bump size for each coefficient
      for (int c=0; c<ncoeff; ++c) {
        if ( naccept[c] > (0.6 * NSTEP) )
          v[c] *= (1 + STEPVAR * ((double)naccept[c] / NSTEP - 0.6) / 0.4);
        else if ( naccept[c] < (0.4 * NSTEP) )
          v[c] /= (1 + STEPVAR * (0.4 - (double)naccept[c] / NSTEP) / 0.4);
        naccept[c] = 0;
      }

      print_simann(k, T, m+1, F, Fopt, group, pot_idx);

      // Rescale potential if necessary
      if ((m+1)%10 == 0 && rescale_me == 1 && mmz->potlist->get_rescale() != 0) {
        Potential *trial_pot = pot->clone();
        if (trial_pot->rescale(comm_inter, out_, 0)) {
          double F2 = trial_pot->compute(comm_inter);
          if (out_) {
            *out_ << "G#" << group << " P#" << pot_idx << " - ";
            *out_ << "F before rescale = " << std::fixed << F << std::endl;
            *out_ << "G#" << group << " P#" << pot_idx << " - ";
            *out_ << "F after rescale  = " << std::fixed << F2 << std::endl;
          }
          int accept = metropolis(pot, trial_pot, &F, F2, T, opt_pot, &Fopt);
          if (accept == 1 && out_) {
            *out_ << "G#" << group << " P#" << pot_idx << " - ";
            *out_ << "Accepted rescale." << std::endl;
          } else if (comm_inter.is_root() && out_) {
            *out_ << "G#" << group << " P#" << pot_idx << " - ";
            *out_ << "Rejected rescale." << std::endl;
          }
        }
        delete trial_pot; // clean up
      }

    } // END STEPS AT FIXED TEMPERATURE

    timer_simann.stop();  // stop the timer
    if (out_) {
      auto ncalls_end = Potential::get_ncalls();
      double time = timer_simann.get_time();  // get the recorded time from start to stop
      *out_ << "G#" << group << " P#" << pot_idx << " - ";
      *out_ << "TTT - calls = " << ncalls_end - ncalls_begin
            << ", total time = " << std::fixed << time
            << " sec, calls/sec = " << (ncalls_end - ncalls_begin)*1.0/time << std::endl;
    }


    T *= TEMPVAR;   // Adjust temperature
    ++k;            // Increment temperature counter
    loopagain = 0;  // Drop out of loop if none of following conditions satisfied

    Fvar[k + NEPS] = F; // record F in backlog of F's
    for (int i = 1; i <= NEPS; ++i)
      if (std::abs(F - Fvar[k + NEPS - i]) > (EPS * F * 0.01)) loopagain = 1;

    // Make sure we don't get stuck in local min
    if (!loopagain && ((F-Fopt) > (EPS*F*0.01))) {
      delete pot;
      pot = opt_pot->clone(); // set back to supposed global min
      F = Fopt;
      loopagain = 1;
      rescale_me = 0;  // turn off rescaling
    }

  } while (k<max_steps && loopagain);

  delete pot; // clean up

  return opt_pot;
}

/* ----------------------------------------------------------------------
   print output from simulated annealing
------------------------------------------------------------------------- */

void OptSimAnn::print_simann(int k, double T, int m, double F, double Fopt, int group, int pot_idx) const
{
  if (!mmz->universe->comm_inter.is_root() || !out_) return;
  if ( k == -1 ) {
    *out_ << std::setw(10) << " ";
    *out_ << std::setw(3) << std::right << "k" << std::setw(5) << " ";
    *out_ << std::setw(12) << std::left << "T";
    *out_ << std::setw(5) << std::right << "m" << std::setw(4) << " ";
    *out_ << std::setw(20) << std::left << "F" << std::setw(3) << " ";
    *out_ << std::setw(1) << std::left << "Fopt";
  } else {
    *out_ << "G#" << group << " P#" << pot_idx << " - ";
    *out_ << std::setw(3) << std::right << k << std::setw(5) << " ";
    *out_ << std::setw(12) << std::left << std::fixed << T;
    *out_ << std::setw(5) << std::right << m << std::setw(4) << " ";
    *out_ << std::setw(20) << std::left << std::scientific << F << std::setw(3) << " ";
    *out_ << std::setw(1) << std::left << std::scientific << Fopt;
  }
  *out_ << std::endl;

  return;
}

/* ----------------------------------------------------------------------
   do a metropolis step: 1 if accepted, 0 if rejected
------------------------------------------------------------------------- */

int OptSimAnn::metropolis(Potential*& pot, Potential *trial_pot, double *F, double F2,
                          double T, Potential*& opt_pot, double *Fopt) const
{
  if ( F2 <= *F ) { // accept this potential
    delete pot;
    pot = trial_pot->clone();
    *F = F2;

    if ( F2 <= *Fopt ) { // record new optimal potential
      delete opt_pot;
      opt_pot = trial_pot->clone();
      *Fopt = F2;
      opt_pot->write_pot(mmz->universe->comm_inter, mmz->potlist->get_tmp_file(), opt_pot->get_global_idx());
    }

    return 1;

  } else if (mmz->random->unif_dist() < std::exp((*F - F2)/T)) { // accept this potential
    delete pot;
    pot = trial_pot->clone();
    *F = F2;
    return 1;
  }

  return 0;
}
