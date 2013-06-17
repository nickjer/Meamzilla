#ifndef MEAMZ_POTENTIAL_H
#define MEAMZ_POTENTIAL_H

#include "pot_fns.h"
#include "pointers.h"
#include "comm.h"

#include <ostream>
#include <vector>

namespace MEAMZ_NS
{

typedef std::string String;
typedef std::vector<String> StringList;

typedef std::vector<double> ErrorVec;                 // Class defining a vector of |error_i|^2

template<typename T>
using Vector = std::vector<T>;

typedef unsigned long long int BigNumber;

class Potential : protected Pointers
{
public:
  enum class Output;                          // Class defining the outputs for the potential

  Potential(class Meamzilla *, int);
  virtual ~Potential();

  virtual Potential* clone() const = 0;       // Copy polymorphic objects

  virtual void init();                        // Initialize potential
  int dbl_cnt_pair() const;                   // Double count pairs?
  double get_max_rcut() const;                // Get maximum radial cutoff from all potentials
  void set_global_idx(int);                   // Set idx position in pot list
  int get_global_idx() const;                 // Get idx position in pot list
  double get_error() const;                   // Get error_sum of potential
  static void reset_ncalls();                 // Reset the number of the times the compute function was called
  static BigNumber get_ncalls();              // Get the number of times the compute function was called
  static void set_output_flag(const Output);  // Set output flag for output
  PotFns& at(int);                            // Get potential function at index value

  int read_pot(const StringList&, int);       // Read in a potential from a potential file
  std::ostream& write(std::ostream&) const;   // Write out object to stream

  int get_ncoeff() const;                     // Return number of adjustable coefficients for this pot
  double& coeff(int);                         // Return value of specified coeffient
  void modify_coeff(int, double, double);     // Perform modification to potential at specified coefficient

  virtual void check_if_similar(Potential *) const = 0; // Check that a potential is similar enough to this
                                                        // potential for pairs/triplets to have identical info

  virtual void compute_trap(const Comm&, int);                          // Trap non-root procs in infinite loop
  virtual double compute(const Comm&, ErrorVec * = nullptr) = 0;        // Compute error sum/vector
  virtual int rescale(const Comm&, std::ostream *, int = 0);            // Rescale potential

  void set_dir(const Vector<double>&, double);                          // Set potential's coeff's to coefficients at user specified direction
  double compute_dir(const Comm&, ErrorVec *, Vector<double>&, double); // Compute error sum at a certain distance in a specific direction from this pot

  void write_pot(const Comm&, String, int = -1) const;    // Write out a single potential to file
  void write_data(const Comm&, String, int = -1);         // Write out fitting data
  void communicate(const Comm&, int);                     // Communicate potential data to other procs in group

  bool operator<(const Potential&) const;     // Compare two potentials using error_sum as comparison

protected:
  int ncoeff_;                        // Number of coefficients to be optimized
  int ntypes_;                        // Number of atom types that this potential describes
  int dbl_cnt_pair_;                  // Whether to double count pairs or not for this potential
  double error_sum_;                  // Error sum of this potential = sum(|error_i|^2,i=0..N-1)
  int global_idx_;                    // Index location of this potential in the global potential list
  static BigNumber ncalls_;           // Number of times the compute fn is called

  int npot_fns_;                      // Number of potential fns in potential (e.g., EAM = 3)
  std::vector<PotFns> pot_fns_;       // List of potential fns (e.g., EAM = phi, rho, U)

  static int is_trapped_;             // Whether procs working on this potential are trapped in infinite loop

  void initialize_pot(const Comm&, ErrorVec * = nullptr);       // Initialize potential on all procs
  void initialize_compute(const Comm&);                         // Initialize computation by setting forces to zero
  void accumulate_error(const Comm&, ErrorVec *, ErrorVec&);    // Construct an error vector and sum up its values squared

  void write_energy(const ErrorVec&, String, int) const;  // Write out energies
  void write_stress(const ErrorVec&, String, int) const;  // Write out stresses
  void write_forces(const ErrorVec&, String, int) const;  // Write out forces
  void write_errors(const ErrorVec&, String, int) const;  // Write out accumulated errors

  virtual void write_punish(const ErrorVec&, String, int) const;  // Write punishment/constraint errors
  virtual void write_extras(const Comm&, String, int);            // Write out extra data

  virtual std::ostream& write_pots(std::ostream&) const;  // Write out actual potential
  virtual std::ostream& write_lmps(std::ostream&) const;  // Write out potential in LAMMPS format

private:
  static Output output_flag_;         // Output flag for output
};

std::ostream& operator<<(std::ostream&, const Potential::Output);
std::ostream& operator<<(std::ostream&, const Potential&);

// Define the Output's of output for the potential
enum class Potential::Output
{
  Pots,   // actual potential
  Lmps    // lammps output of actual potential
};

} // namespace MEAMZ_NS

#endif // MEAMZ_POTENTIAL_H
