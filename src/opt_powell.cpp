/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "opt_powell.h"
#include "input.h"
#include "random.h"
#include "timer.h"
#include "error.h"
#include "line_min_brent.h"

#include <cmath>
#include <iomanip>
#include <functional>

extern "C" {
void dsysvx_(char *fact, char *uplo, int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf,
             int *ipiv, double *b, int *ldb, double *x, int *ldx, double *rcond, double *ferr,
             double *berr, double *work, int *lwork, int *iwork, int *info);
}

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

OptPowell::OptPowell(class Meamzilla *mmz, std::ostream *out) : Optimizer(mmz, out), d_eps_(0.0)
{
  max_steps_ = 1000000;
}

/* ---------------------------------------------------------------------- */

OptPowell::~OptPowell()
{
  //dtor
}

/* ----------------------------------------------------------------------
   initialize optimizer
------------------------------------------------------------------------- */

void OptPowell::init()
{
  if (mmz->input) {
    mmz->input->parse("d_eps", 0, d_eps_);
    mmz->input->parse("max_steps", 0, max_steps_);
  }

  return;
}

/* ----------------------------------------------------------------------
   use the optimizer on a potential list
------------------------------------------------------------------------- */

void OptPowell::compute(PotList& pot_list)
{
  if (mmz->universe->is_root() && out_) {
    *out_ << "Beginning Powell least squares minimization..." << std::endl << std::endl;
    print_powell(-1, 0, 0, 0, 0);
  }
  mmz->universe->barrier();

  // Loop through all the potentials this group sees
  int npots = pot_list.get_group_npots();
  for (int p=0; p<npots; ++p) {
    mmz->random->reset_seed();    // reset the seed to original seed for reproducibility
    pot_list[p]->compute_trap(mmz->universe->comm_inter, 1);  // trap procs that aren't root of comm_inter
    if (!mmz->universe->comm_inter.is_root()) break;
    Potential *pot = compute_pot(pot_list[p], max_steps_);
    pot_list[p]->compute_trap(mmz->universe->comm_inter, 0);  // untrap procs
    delete pot_list[p];
    pot_list[p] = pot;
  } // End loop over potentials

  pot_list.write_pots(pot_list.get_end_file() + ".powell"); // this is also a barrier for procs and passes all pots to all root nodes

  if (mmz->universe->is_root() && out_)
    *out_ << std::endl << "Finished Powell least squares minimization..." << std::endl << std::endl;
  mmz->universe->barrier();

  return;
}

/* ----------------------------------------------------------------------
   optimization loop
------------------------------------------------------------------------- */

Potential* OptPowell::compute_pot(const Potential *main_pot, int max_steps) const
{
  // Some local variables to make my life easier
  Comm& comm_inter = mmz->universe->comm_inter;     // Communicator for procs in a group
  int group = mmz->universe->get_group();           // Group number
  int pot_idx = main_pot->get_global_idx();         // Global potential index in list of potentials for this pot
  int ncoeff = main_pot->get_ncoeff();              // Number of coefficients to be optimized
  ErrorVec error_vec;                               // Vector of fitting errors

  Potential *pot = main_pot->clone();               // Clone main potential to local potential
  double F = pot->compute(comm_inter, &error_vec);  // Get force error
  int esize = error_vec.size();                     // Number of functions (fitting errors) in least squares method

  if (max_steps == 0) return pot;

  // Initialize matrices
  Matrix<double> gamma(esize, ncoeff);              // Matrix of derivatives (esize x ncoeff)
  Matrix<double> d(ncoeff, ncoeff);                 // Direction vectors (ncoeff x ncoeff)
  Matrix<double> lineqsys(ncoeff, ncoeff);          // Lin.Eq.Sys. Matrix (ncoeff x ncoeff)
  Vector<double> p(ncoeff);                         // Vectors needed in Powell's algorithm (ncoeff)
  Vector<double> q(ncoeff);
  Vector<double> delta(ncoeff);                     // Vector pointing into correct direction (ncoeff)
  double delta_mag;                                 // Magnitude of delta vector

  double F_older = F;                               // used for storing outer loop F
  double dF;                                        // used for difference between F's in inner loop

  if (F < NOTHING) {   // Already found minimum F
    if (out_)
      *out_ << "Error already too small to optimize, aborting..." << std::endl;
    return pot;
  }

  // Begin outer loop
  Timer timer_powell;  // setup timer and start it
  auto ncalls_outer_begin = Potential::get_ncalls();  // Get number of compute calls for potential

  int n = 0;
  do {
    auto ncalls_inner_begin = Potential::get_ncalls();  // Get number of compute calls for potential

    // Initialize the gamma matrix and direction "d" vectors
    if (int i = gamma_init(pot, error_vec, esize, ncoeff, &d, &gamma)) {
      if (mmz->potlist->get_rescale() != 0 && out_)
        *out_ << "F does not depend on coefficient " << i-1 <<
                 ", trying to rescale if it is an option!" << std::endl;

      if (mmz->potlist->get_rescale() != 0) {
        pot->rescale(comm_inter, out_, 1);
        F = pot->compute(comm_inter, &error_vec);
      }

      if (gamma_init(pot, error_vec, esize, ncoeff, &d, &gamma)) {
        pot->write_pot(comm_inter, mmz->potlist->get_tmp_file(), pot_idx);
        if (out_)
          *out_ << "F does not depend on coefficient " << i-1 << ", fit impossible!" << std::endl;
        break;
      }
    }

    // Initialize the linear equation system matrices (lineqsys and p)
    lineqsys_init(error_vec, esize, ncoeff, gamma, &lineqsys, &p);

    F_older = F;  // store outer loop F

    // Begin inner loop
    int m = 0;
    do {
      double F_old = F; // store inner loop F

      // (a) Solve linear equations

      double cond;             // Condition number dsysvx
      double ferror, berror;   // forward/backward error estimates
      int info;                // Info integer
      char uplo = 'U';         // Upper triangle of A is stored
      char fact = 'N';         // Specifies whether or not the factored form of A has been supplied on entry.
      int nrhs = 1;            // The number of right hand sides, i.e., the number of columns of the matrix B
      Vector<int> perm_indx(ncoeff);              // Keeps track of LU pivoting
      Matrix<double> les_inverse(ncoeff, ncoeff); // LU decomp. of the lineqsys
      int worksize = 64 * ncoeff;
      Vector<double> work(worksize);
      Vector<int> iwork(ncoeff);

      dsysvx_(&fact, &uplo, &ncoeff, &nrhs, &lineqsys(0,0), &ncoeff, &les_inverse(0,0), &ncoeff,
              &perm_indx[0], &p[0], &ncoeff, &q[0], &ncoeff, &cond, &ferror, &berror,
              &work[0], &worksize, &iwork[0], &info);

      if (info > 0 && info <= ncoeff) {
        if (out_)
          *out_ << "Linear equation system singular after step " << m << " i=" << info << std::endl;
        break;
      }

      // (b) Get delta by multiplying q with the direction vectors

      delta_mag = 0.0;
      for (int i=0; i<ncoeff; ++i) {
        delta[i] = 0.;
        for (int j=0; j<ncoeff; ++j) {
          delta[i] += d(i,j) * q[j];
        }
        delta_mag += delta[i]*delta[i];
      }

      // (c) Minimize F(coeff's) along vector delta, return new F
      // * Modifies potential pointer passed to it

      LineMinBrent linmin;
      try {
        std::function<double(double)> eval = std::bind(&Potential::compute_dir, pot,
                                                       std::cref(comm_inter), nullptr, std::ref(delta), std::placeholders::_1);
        linmin.set_func(eval);
        F = linmin.find_minima();
        double x_min = linmin.get_xmin1();
        pot->set_dir(delta, x_min);
        pot->compute(comm_inter, &error_vec);

      } catch(Error& ex) {
        if (out_)
          *out_ << ex.what() << "pot=" << pot_idx << std::endl;
        mmz->universe->comm_all.abort_one();
      }

      // (d) If error estimate is too high after minimization in 5 directions: restart outer loop

      if (ferror+berror>1. && m>5) break;

      // (e) Find optimal direction to replace: p_i * q_i =(approx) delta(F) in that direction
      // so it is the component vector with largest change in F that will be replaced in direction vectors

      int max_col = 0;
      double max_change = 0.;
      for (int i=0; i<ncoeff; ++i) {
        if (std::abs(p[i]*q[i]) > max_change) {
          max_change = std::abs(p[i]*q[i]);
          max_col = i;
        }
      }

      // (f) update gamma, but if fn returns 1, matrix will be singular, break inner loop and restart with new matrix
      // slao normalize delta in this update so it is properly scaled for part (g)
      // Modifies: gamma and delta only

      if (gamma_update(max_col, pot, error_vec, esize, ncoeff, linmin, &delta, &gamma)) {
        if (out_)
          *out_ << "Matrix gamma singular after step " << m << ", restarting inner loop" << std::endl;
        break;
      }

      // (g) Set new *scaled* direction vector

      for (int i=0; i<ncoeff; ++i)
        d(i,max_col) = delta[i];

      // (h) update linear equation system
      // Modifies: lineqsys and p only

      lineqsys_update(max_col, error_vec, esize, ncoeff, gamma, &lineqsys, &p);

      dF = F_old - F;

      ++m;  // increment inner loop
    } while ((m<ncoeff || (m<=INNERLOOPS && dF>PRECISION)) && dF<TOOBIG && delta_mag<TOOBIG);  // END INNER LOOP

    print_powell(m, F, Potential::get_ncalls() - ncalls_inner_begin, group, pot_idx);

    // Rescale every 4th step
    if ((n+1)%4 == 0 && mmz->potlist->get_rescale() != 0) {
      Potential *trial_pot = pot->clone();
      if (trial_pot->rescale(comm_inter, out_, 0)) {
        double F2 = trial_pot->compute(comm_inter);
        if (out_) {
          *out_ << "G#" << group << " P#" << pot_idx << " - ";
          *out_ << "F before rescale = " << std::fixed << F << std::endl;
          *out_ << "G#" << group << " P#" << pot_idx << " - ";
          *out_ << "F after rescale  = " << std::fixed << F2 << std::endl;
        }
        int accept = 0;
        if (F2 <= F) {  // accept potential
          // swap potentials
          delete pot;
          pot = trial_pot;
          trial_pot = nullptr;
          F = pot->compute(comm_inter, &error_vec);
          accept = 1;
        }
        if (accept == 1 && out_) {
          *out_ << "G#" << group << " P#" << pot_idx << " - ";
          *out_ << "Accepted rescale." << std::endl;
        } else if (out_) {
          *out_ << "G#" << group << " P#" << pot_idx << " - ";
          *out_ << "Rejected rescale." << std::endl;
        }
        if (trial_pot) delete trial_pot;
      }
    }

    if ((n+1)%10 == 0) {
      timer_powell.stop();  // stop the timer
      auto ncalls_end = Potential::get_ncalls();
      double time = timer_powell.get_time();  // get the recorded time from start to stop
      if (out_) {
        *out_ << "G#" << group << " P#" << pot_idx << " - ";
        *out_ << "TTT - calls = " << ncalls_end - ncalls_outer_begin
              << ", total time = " << std::fixed << time
              << " sec, calls/sec = " << (ncalls_end - ncalls_outer_begin)*1.0/time << std::endl;
      }
      timer_powell.start();
      ncalls_outer_begin = Potential::get_ncalls();
    }

    pot->write_pot(comm_inter, mmz->potlist->get_tmp_file(), pot_idx);

    ++n;  // increment outer loop
  } while (((F_older - F > PRECISION/10.) || (F_older - F < 0))
           && (F_older - F > d_eps_) && n < max_steps);  // END OUTER LOOP

  return pot;
}

/* ----------------------------------------------------------------------
   print output from simulated annealing
------------------------------------------------------------------------- */

void OptPowell::print_powell(int m, double F, int ncalls, int group, int pot_idx) const
{
  if (!mmz->universe->comm_inter.is_root() || !out_) return;

  if ( m == -1 ) {
    *out_ << std::setw(10) << " ";
    *out_ << std::setw(3) << std::right << "n";
    *out_ << std::setw(8) << std::right << "ncalls" << "   ";
    *out_ << "F";
  } else {
    *out_ << "G#" << group << " P#" << pot_idx << " - ";
    *out_ << std::setw(3) << std::right << m;
    *out_ << std::setw(8) << std::right << ncalls << "   ";
    *out_ << std::scientific << std::setprecision(8) << F;
  }
  *out_ << std::endl;

  return;
}

/* ----------------------------------------------------------------------
   gamma_init: (Re-)Initialize gamma[j][i] (Gradient Matrix) after
   (Re-)Start or whenever necessary by calculating numerical
   gradients in coordinate directions. Includes re-setting the
   direction vectors to coordinate directions.
------------------------------------------------------------------------- */

int OptPowell::gamma_init(const Potential *pot, const ErrorVec& evec, int esize,
                          int ncoeff, Matrix<double> *d_ptr, Matrix<double> *gamma_ptr) const
{
  // Tmp variables to make life easier
  Potential *tmp_pot = pot->clone();
  Matrix<double>& d = *d_ptr;
  Matrix<double>& gamma = *gamma_ptr;

  // Initialize direction vectors to coordinate directions
  // d_ij = KroneckerDelta_ij
  for (int i=0; i<ncoeff; ++i)
    for (int j=0; j<ncoeff; ++j)
      d(i,j) = (i == j) ? 1. : 0.;

  // Initialize gamma by calculating numerical derivatives
  for (int i=0; i<ncoeff; ++i) {
    ErrorVec tmp_evec;
    double store_coeff = tmp_pot->coeff(i);
    tmp_pot->coeff(i) += EPS; // increase coefficient by small amount
    tmp_pot->compute(mmz->universe->comm_inter, &tmp_evec);

    double sum = 0.0;
    for (int j=0; j<esize; ++j) {
      double temp = (tmp_evec[j] - evec[j])/EPS;
      gamma(j,i) = temp;
      sum += temp * temp;
    }

    double mag_i = std::sqrt(sum);

    tmp_pot->coeff(i) = store_coeff;  // reset coefficient

    if (mag_i < NOTHING) return i+1;  // this is singular matrix so abort

    for (int j=0; j<esize; ++j) gamma(j,i) /= mag_i; // normalize this gamma column
    d(i,i) /= mag_i; // rescale direction vector as well
  }

  delete tmp_pot;

  return 0;
}

/* ----------------------------------------------------------------------
   gamma_update: Update column max_col of gamma (to newly calculated
   numerical derivatives (calculated from f1, f2 at lambda1,lambda2);
   normalize new vector
------------------------------------------------------------------------- */

int OptPowell::gamma_update(int max_col, Potential *pot, const ErrorVec& evec, int esize, int ncoeff,
                            LineMinBrent& linmin, Vector<double> *delta_ptr, Matrix<double> *gamma_ptr) const
{
  double mu = 0.0;
  Vector<double>& delta = *delta_ptr;
  Matrix<double>& gamma = *gamma_ptr;

  double lambda1 = linmin.get_xmin1();
  double lambda2 = linmin.get_xmin2();
  double F = linmin.get_Fmin1();
  const ErrorVec& f1 = evec;
  ErrorVec f2;
  pot->compute_dir(mmz->universe->comm_inter, &f2, delta, lambda2-lambda1); // pot is already displaced by lambda1 in direction

  for (int j=0; j<esize; ++j) {
    double temp = (f1[j] - f2[j])/(lambda1 - lambda2);
    gamma(j,max_col) = temp;
    mu += temp * f1[j];
  }
  mu /= F;

  // Similar to gamma_init now
  double sum = 0.0;
  for (int j=0; j<esize; ++j) {
    double temp = gamma(j,max_col)- mu * f1[j];
    gamma(j,max_col) = temp;
    sum += temp * temp;
  }

  double mag = std::sqrt(sum);

  if (mag < NOTHING) return 1;

  for (int j=0; j<esize; ++j) gamma(j,max_col) /= mag;
  for (int i=0; i<ncoeff; ++i) delta[i] /= mag;

  return 0;
}

/* ----------------------------------------------------------------------
   lineqsys_init: Initialize LinEqSys matrix, vector p in
                      lineqsys . q == p
------------------------------------------------------------------------- */

void OptPowell::lineqsys_init(const ErrorVec& evec, int esize, int ncoeff, const Matrix<double>& gamma,
                              Matrix<double> *lineqsys_ptr, Vector<double> *p_ptr) const
{
  // Tmp variables to make life easier
  Matrix<double>& lineqsys = *lineqsys_ptr;
  Vector<double>& p = *p_ptr;

  // Compute p vector
  for (int i=0; i<ncoeff; ++i) {
    p[i] = 0.;
    for (int j=0; j<esize; ++j) p[i] -= gamma(j,i) * evec[j];
  }

  // Compute lineqsys matrix (starting with diagonal components)
  for (int i=0; i<ncoeff; ++i) {
    lineqsys(i,i) = 0.;
    for (int j=0; j<esize; ++j) lineqsys(i,i) += gamma(j,i) * gamma(j,i);

    for (int j=i+1; j<ncoeff; ++j) {
      lineqsys(i,j) = 0.;
      for (int k=0; k<esize; ++k) lineqsys(i,j) += gamma(k,i) * gamma(k,j);
      lineqsys(j,i) = lineqsys(i,j);
    }
  }

  return;
}

/* ----------------------------------------------------------------------
   lineqsys_update: Update LinEqSys matrix row and
   column max_col in vector p.
------------------------------------------------------------------------- */

void OptPowell::lineqsys_update(int max_col, const ErrorVec& evec, double esize, double ncoeff,
                                const Matrix<double>& gamma, Matrix<double> *lineqsys_ptr, Vector<double> *p_ptr) const
{
  // Tmp variables to make life easier
  Matrix<double>& lineqsys = *lineqsys_ptr;
  Vector<double>& p = *p_ptr;

  for (int i=0; i<ncoeff; ++i) {
    p[i] = 0.;
    for (int j=0; j<esize; ++j) p[i] -= gamma(j,i) * evec[j];
  }

  for (int j=0; j<ncoeff; ++j) {
    lineqsys(max_col,j) = 0.;
    for (int k=0; k<esize; ++k) lineqsys(max_col,j) += gamma(k,max_col) * gamma(k,j);
    lineqsys(j,max_col) = lineqsys(max_col,j);
  }

  return;
}
