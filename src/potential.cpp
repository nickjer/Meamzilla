#include "potential.h"
#include "style_basis.h"
#include "universe.h"
#include "config.h"
#include "cell.h"
#include "atom.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace MEAMZ_NS;

// Initialize potential output flag
Potential::Output Potential::output_flag_ = Potential::Output::Pots;

// Initialize number of compute calls to 0
BigNumber Potential::ncalls_ = 0;

// Initialize whether procs are trapped
int Potential::is_trapped_ = 0;

/* ---------------------------------------------------------------------- */

Potential::Potential(class Meamzilla *mmz, int ntypes) : Pointers(mmz), ncoeff_(0),
                                                         dbl_cnt_pair_(0), error_sum_(0),
                                                         npot_fns_(0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

Potential::~Potential()
{
  //dtor
}

/* ----------------------------------------------------------------------
   initialize potential
------------------------------------------------------------------------- */

void Potential::init()
{
  return;
}

/* ----------------------------------------------------------------------
   double count pairs?
------------------------------------------------------------------------- */

int Potential::dbl_cnt_pair() const
{
  return dbl_cnt_pair_;
}

/* ----------------------------------------------------------------------
   get maximum radial cutoff from all potentials
------------------------------------------------------------------------- */

double Potential::get_max_rcut() const
{
  double max_rcut = 0.0;

  for (PotFns pot_fn : pot_fns_)
    if (pot_fn.is_radial())
      max_rcut = std::max(max_rcut, pot_fn.get_max_rcut());

  return max_rcut;
}

/* ----------------------------------------------------------------------
   set idx position in pot list
------------------------------------------------------------------------- */

void Potential::set_global_idx(int idx)
{
  global_idx_ = idx;
  return;
}

/* ----------------------------------------------------------------------
   get idx position in pot list
------------------------------------------------------------------------- */

int Potential::get_global_idx() const
{
  return global_idx_;
}

/* ----------------------------------------------------------------------
   get error_sum of potential
------------------------------------------------------------------------- */

double Potential::get_error() const
{
  return error_sum_;
}

/* ----------------------------------------------------------------------
   reset the number of the times the compute function was called
------------------------------------------------------------------------- */

void Potential::reset_ncalls()
{
  ncalls_ = 0;
}

/* ----------------------------------------------------------------------
   get the number of times the compute function was called
------------------------------------------------------------------------- */

unsigned long long int Potential::get_ncalls()
{
  return ncalls_;
}

/* ----------------------------------------------------------------------
   set output flag for output
------------------------------------------------------------------------- */

void Potential::set_output_flag(const Output flag)
{
  output_flag_ = flag;
  return;
}

/* ----------------------------------------------------------------------
   get potential function at index value
------------------------------------------------------------------------- */

PotFns& Potential::at(int idx)
{
  return pot_fns_[idx];
}

/* ----------------------------------------------------------------------
   read in a potential from a potential file
------------------------------------------------------------------------- */

int Potential::read_pot(const StringList& lines, int line_num)
{
  for (PotFns& pot_fn : pot_fns_) {
    line_num = pot_fn.read_pot_fns(lines, line_num);
    ncoeff_ += pot_fn.get_ncoeff();
  }

  return line_num;
}

/* ----------------------------------------------------------------------
   write out to stream
------------------------------------------------------------------------- */

std::ostream& Potential::write(std::ostream& stream) const
{
  switch (output_flag_) {
  case Output::Pots:
    return write_pots(stream);
  case Output::Lmps:
    return write_lmps(stream);
  default:
    return stream;
  }
}

/* ----------------------------------------------------------------------
   return number of adjustable coefficients for this pot
------------------------------------------------------------------------- */

int Potential::get_ncoeff() const
{
  return ncoeff_;
}

/* ----------------------------------------------------------------------
   return reference to coefficient in potential that can be modified
------------------------------------------------------------------------- */

double& Potential::coeff(int idx)
{
  double *value = nullptr;
  for (PotFns& pot_fn : pot_fns_) {
    value = &pot_fn.coeff(idx);  // idx is allowed to change
    if (idx < 0) return *value;
  }

  // If it fails output the error
  std::cerr << "ERROR: unable to access coefficient for optimization" << std::endl;
  mmz->universe->comm_all.abort_one();

  // It will never get this far (hopefully!)
  return *value;
}

/* ----------------------------------------------------------------------
   perform modification to potential at specified coefficient
------------------------------------------------------------------------- */

void Potential::modify_coeff(int idx, double height, double width)
{
  // Find potential function where idx lies in, and modify that pot fn
  int sum = 0;
  for (PotFns& pot_fn : pot_fns_) {
    int pot_fn_ncoeff = pot_fn.get_ncoeff();
    sum += pot_fn_ncoeff;
    if (idx < sum) {
      pot_fn.modify_coeff(idx + pot_fn_ncoeff - sum, height, width);
      return;
    }
  }

  return;
}

/* ----------------------------------------------------------------------
   trap non-root procs in infinite loop
------------------------------------------------------------------------- */

void Potential::compute_trap(const Comm& comm, int flag)
{
  Vector<double> tmp_vec; // nothing is ever written to this variable

  if (comm.is_root() && flag)
    is_trapped_ = 1;
  else if (comm.is_root())
    is_trapped_ = 0;

  do {
    comm.bcast(&flag, 1, MPI_INT, comm.get_root());

    if (flag == 2)
      compute(comm);
    else if (flag == 3)
      compute(comm, &tmp_vec);

  } while (flag && !comm.is_root());

  return;
}

/* ----------------------------------------------------------------------
   rescale potential
------------------------------------------------------------------------- */

int Potential::rescale(const Comm& comm, std::ostream *out, int flag)
{
  return 0;
}

/* ----------------------------------------------------------------------
   set potential's coeff's to coefficients at user
   specified direction
------------------------------------------------------------------------- */

void Potential::set_dir(const Vector<double>& dir_vec, double x)
{
  for (int i=0; i<ncoeff_; ++i)
    coeff(i) += x * dir_vec[i];

  return;
}

/* ----------------------------------------------------------------------
   compute error sum at a certain distance in a specific
   direction from this pot
------------------------------------------------------------------------- */

double Potential::compute_dir(const Comm& comm, ErrorVec *err_vec, Vector<double>& dir_vec, double x)
{
  Potential *tmp_pot = clone();

  tmp_pot->set_dir(dir_vec, x);

  double err_sum = tmp_pot->compute(comm, err_vec);

  delete tmp_pot;

  return err_sum;
}


/* ----------------------------------------------------------------------
   write out a single potential to file
------------------------------------------------------------------------- */

void Potential::write_pot(const Comm& comm, String filename, int num) const
{
  if (!comm.is_root()) return;

  std::ostringstream oss;
  oss << filename;
  if (num >= 0)
    oss << "." << num;

  filename = oss.str();

  std::ofstream ofs;
  ofs.open(filename.c_str());

  if (ofs)
    ofs << Output::Pots << *this;

  ofs.close();

  return;
}

/* ----------------------------------------------------------------------
   write out fitting data
------------------------------------------------------------------------- */

void Potential::write_data(const Comm& comm, String filename, int num)
{
  ErrorVec e_vec;
  compute(comm, &e_vec);

  write_extras(comm, filename, num);

  if (comm.is_root()) {
    write_energy(e_vec, filename, num);
    write_stress(e_vec, filename, num);
    write_forces(e_vec, filename, num);
    write_errors(e_vec, filename, num);

    write_punish(e_vec, filename, num);
  }

  return;
}

/* ----------------------------------------------------------------------
   communicate potential data to other procs in group
------------------------------------------------------------------------- */

void Potential::communicate(const Comm& comm, int root)
{
  comm.bcast(&error_sum_, 1, MPI_DOUBLE, root);

  for (PotFns& pot_fn : pot_fns_)
    pot_fn.communicate(comm, root);

  return;
}

/* ----------------------------------------------------------------------
   compare two potentials using error_sum as comparison
------------------------------------------------------------------------- */

bool Potential::operator<(const Potential& rhs) const
{
  return error_sum_ < rhs.error_sum_;
}

/* ----------------------------------------------------------------------
   construct an error vector and sum up its values squared
------------------------------------------------------------------------- */

void Potential::initialize_pot(const Comm& comm, ErrorVec *error_vec)
{
  // Reset 2nd derivs for all potentials
  if (comm.is_root())
    for (PotFns& pot_fn : pot_fns_)
      pot_fn.refresh_basis();

  // Broadcast this potential to procs
  communicate(comm, comm.get_root());

  return;
}

/* ----------------------------------------------------------------------
   initialize potential by resetting forces
------------------------------------------------------------------------- */

void Potential::initialize_compute(const Comm& comm)
{
  // Reset 2nd derivs for all potentials
  if (comm.is_root())
    for (PotFns& pot_fn : pot_fns_)
      pot_fn.refresh_basis();

  // Broadcast this potential to procs
  communicate(comm, comm.get_root());

  // Reset forces, energies, stresses, and densities for every atom
  for (Cell*& cell_ptr : mmz->config->cells) {
    cell_ptr->energy = 0.0;
    cell_ptr->stress = 0.0;

    for (Atom*& atom_ptr : cell_ptr->atoms)
      atom_ptr->force = 0.0;
  }

  return;
}


/* ----------------------------------------------------------------------
   construct an error vector and sum up its values squared
------------------------------------------------------------------------- */

void Potential::accumulate_error(const Comm& comm, ErrorVec *error_vec, ErrorVec& addon_error)
{
  ///////////////////////////////////////////////////
  // PACK PROPERTIES OVER MPI
  ///////////////////////////////////////////////////

  // Define temporary error vector
  // (1 energy, 6 stresses)*ncells = 7*ncells
  // (3 forces per atom)*natoms = 3*natoms
  int tmp_err_size = 7*mmz->config->ncells + 3*mmz->config->total_natoms;
  ErrorVec tmp_err_vec(tmp_err_size);

  // Addon extra user-supplied error
  tmp_err_size += addon_error.size();
  tmp_err_vec.insert(tmp_err_vec.end(), addon_error.begin(), addon_error.end());

  // Fill in the error vector
  int counter = -1;
  for (Cell*& cell_ptr : mmz->config->cells) {
    tmp_err_vec[++counter]   = cell_ptr->energy;
    for (int j=0; j<6; ++j)
      tmp_err_vec[++counter] = cell_ptr->stress[j];
    for (Atom*& atom_ptr : cell_ptr->atoms) {
      tmp_err_vec[++counter] = atom_ptr->force.x;
      tmp_err_vec[++counter] = atom_ptr->force.y;
      tmp_err_vec[++counter] = atom_ptr->force.z;
    }
  }

  // Sum up these values across processors
  ErrorVec tmp_error_final(tmp_err_size, 0.0);
  comm.reduce(&tmp_err_vec[0],&tmp_error_final[0],tmp_err_size,MPI_DOUBLE,MPI_SUM,comm.get_root());

  ///////////////////////////////////////////////////
  // UNPACK PROPERTIES OVER MPI
  ///////////////////////////////////////////////////

  error_sum_ = 0.0; // Set error

  counter = -1;
  for (Cell*& cell_ptr : mmz->config->cells) {
    cell_ptr->energy       = tmp_error_final[++counter]/cell_ptr->natoms;
    tmp_err_vec[counter]   = cell_ptr->weight * cell_ptr->eweight * mmz->config->get_energy_weight() *
                             cell_ptr->natoms * (cell_ptr->energy - cell_ptr->energy0);
    error_sum_        += tmp_err_vec[counter]*tmp_err_vec[counter];

    for (int j=0; j<6; ++j) {
      cell_ptr->stress[j]  = tmp_error_final[++counter]/cell_ptr->vol;
      tmp_err_vec[counter] = cell_ptr->weight * cell_ptr->sweight * mmz->config->get_stress_weight() *
                            (cell_ptr->stress[j] - cell_ptr->stress0[j]);
      error_sum_      += tmp_err_vec[counter]*tmp_err_vec[counter];
    }

    for (Atom*& atom_ptr : cell_ptr->atoms) {
      atom_ptr->force.x    = tmp_error_final[++counter];
      tmp_err_vec[counter] = cell_ptr->weight * (atom_ptr->force.x - atom_ptr->force0.x);
      error_sum_      += tmp_err_vec[counter]*tmp_err_vec[counter];
      atom_ptr->force.y    = tmp_error_final[++counter];
      tmp_err_vec[counter] = cell_ptr->weight * (atom_ptr->force.y - atom_ptr->force0.y);
      error_sum_      += tmp_err_vec[counter]*tmp_err_vec[counter];
      atom_ptr->force.z    = tmp_error_final[++counter];
      tmp_err_vec[counter] = cell_ptr->weight * (atom_ptr->force.z - atom_ptr->force0.z);
      error_sum_      += tmp_err_vec[counter]*tmp_err_vec[counter];
    }
  }

  // Unpack addon error now
  while (++counter < tmp_err_size ) {
    tmp_err_vec[counter] = tmp_error_final[counter];
    error_sum_ += tmp_err_vec[counter]*tmp_err_vec[counter];
  }

  if (error_vec && comm.is_root())
    error_vec->swap(tmp_err_vec);

  return;
}

/* ----------------------------------------------------------------------
   write out energy data
------------------------------------------------------------------------- */

void Potential::write_energy(const ErrorVec& e_vec, String filename, int num) const
{
  filename += ".energy";
  std::ostringstream oss;
  oss << filename;
  if (num >= 0)
    oss << "." << num;

  std::ofstream energy_fs (oss.str(), std::ios_base::out);
  energy_fs << "Global energy weight is " << std::fixed << mmz->config->get_energy_weight() << std::endl;
  energy_fs << std::setw(3)  << std::right << "#"
            << std::setw(11) << std::right << "conf_w"
            << std::setw(11) << std::right << "conf_ew"
            << std::setw(9)  << std::right << "natoms"
            << std::setw(15) << std::right << "energy error"
            << std::setw(18) << std::right << "e"
            << std::setw(18) << std::right << "e0"
            << std::setw(15) << std::right << "|e-e0|"
            << std::setw(15) << std::right << "e-e0"
            << std::setw(15) << std::right << "de/e0"
            << std::endl;

  int energy_idx = 0;
  for (Cell*& cell_ptr : mmz->config->cells) {
    Cell& cell = *cell_ptr;
    energy_fs << std::setw(3)  << std::right << cell.cell_idx << std::fixed
              << std::setw(11) << std::right << std::setprecision(4) << cell.weight
              << std::setw(11) << std::right << std::setprecision(4) << cell.eweight
              << std::setw(9)  << std::right << cell.natoms
              << std::setw(15) << std::right << std::setprecision(6) << e_vec[energy_idx]*e_vec[energy_idx]
              << std::setw(18) << std::right << std::setprecision(10) << cell.energy
              << std::setw(18) << std::right << std::setprecision(10) << cell.energy0
              << std::setw(15) << std::right << std::setprecision(6) << std::abs(cell.energy-cell.energy0)
              << std::setw(15) << std::right << std::setprecision(6) << cell.energy-cell.energy0;
    if (cell.energy0==0)
      energy_fs << std::setw(15) << std::right << "inf";
    else
      energy_fs << std::setw(15) << std::right << std::setprecision(6) << (cell.energy0-cell.energy)/cell.energy0;
    energy_fs << std::endl;
    energy_idx += 7 + 3*cell.natoms;
  }

  energy_fs.close();

  return;
}

/* ----------------------------------------------------------------------
   write out stress data
------------------------------------------------------------------------- */

void Potential::write_stress(const ErrorVec& e_vec, String filename, int num) const
{
  filename += ".stress";
  std::ostringstream oss;
  oss << filename;
  if (num >= 0)
    oss << "." << num;

  std::ofstream stress_fs (oss.str(), std::ios_base::out);
  stress_fs << "Global stress weight is " << std::fixed << mmz->config->get_stress_weight() << std::endl;
  stress_fs << std::setw(3)  << std::right << "#"
            << std::setw(11) << std::right << "conf_w"
            << std::setw(11) << std::right << "conf_sw"
            << std::setw(15) << std::right << "stress error"
            << std::setw(18) << std::right << "s"
            << std::setw(18) << std::right << "s0"
            << std::setw(15) << std::right << "ds/s0"
            << std::endl;

  int stress_idx = 0;
  for (Cell*& cell_ptr : mmz->config->cells) {
    Cell& cell = *cell_ptr;
    ++stress_idx;
    for (int dim=0; dim<6; ++dim) {
      stress_fs << std::setw(3)  << std::right << cell.cell_idx << std::fixed
                << std::setw(11) << std::right << std::setprecision(4)  << cell.weight
                << std::setw(11) << std::right << std::setprecision(4)  << cell.sweight
                << std::setw(15) << std::right << std::setprecision(6)  << e_vec[stress_idx]*e_vec[stress_idx]
                << std::setw(18) << std::right << std::setprecision(10) << cell.stress[dim]
                << std::setw(18) << std::right << std::setprecision(10) << cell.stress0[dim];
      if (cell.stress0[dim]==0)
        stress_fs << std::setw(15) << std::right << "inf";
      else
        stress_fs << std::setw(15) << std::right << std::setprecision(6)<< (cell.stress0[dim]-cell.stress[dim])/cell.stress0[dim];
      stress_fs << std::endl;
      ++stress_idx;
    }
    stress_idx += 3*cell.natoms;
  }

  stress_fs.close();

  return;
}

/* ----------------------------------------------------------------------
   write out force data
------------------------------------------------------------------------- */

void Potential::write_forces(const ErrorVec& e_vec, String filename, int num) const
{
  filename += ".force";
  std::ostringstream oss;
  oss << filename;
  if (num >= 0)
    oss << "." << num;

  std::ofstream force_fs (oss.str(), std::ios_base::out);
  force_fs << std::setw(6) << std::right << "conf:"
           << std::setw(5) << std::right << "atom"
           << std::setw(8) << std::right << "type"
           << std::setw(15) << std::right << "force error"
           << std::setw(15) << std::right << "f"
           << std::setw(15) << std::right << "f0"
           << std::setw(15) << std::right << "df/f0"
           << std::endl;

  int force_idx = 0;
  for (Cell*& cell_ptr : mmz->config->cells) {
    Cell& cell = *cell_ptr;
    force_idx += 7; // cycle through 1 energy and 6 stresses per cell
    for (Atom*& atom_ptr : cell.atoms) {
      for (int dim=0; dim<3; ++dim) {
        force_fs << std::setw(5) << std::right << cell.cell_idx << ":"
                 << std::setw(5) << std::right << atom_ptr->atom_idx
                 << std::setw(8) << std::right << atom_ptr->typ
                 << std::setw(15) << std::right << std::fixed << e_vec[force_idx]*e_vec[force_idx]
                 << std::setw(15) << std::right << std::fixed << atom_ptr->force[dim]
                 << std::setw(15) << std::right << std::fixed << atom_ptr->force0[dim];
        if (atom_ptr->force0[dim]==0)
          force_fs << std::setw(15) << std::right << "inf";
        else
          force_fs << std::setw(15) << std::right << std::fixed
                   << (atom_ptr->force0[dim]-atom_ptr->force[dim])/atom_ptr->force0[dim];
        force_fs << std::endl;
        ++force_idx;
      }
    }
  }

  force_fs.close();

  return;
}

/* ----------------------------------------------------------------------
   write out error data
------------------------------------------------------------------------- */

void Potential::write_errors(const ErrorVec& e_vec, String filename, int num) const
{
  filename += ".error";
  std::ostringstream oss;
  oss << filename;
  if (num >= 0)
    oss << "." << num;

  int esize = e_vec.size();
  int ncells = mmz->config->ncells;
  int natoms = mmz->config->total_natoms;

  double sum_energy = 0.0;
  double sum_stress = 0.0;
  double sum_forces = 0.0;
  double sum_punish = 0.0;

  int idx = 0;
  for (Cell*& cell_ptr : mmz->config->cells) {
    Cell& cell = *cell_ptr;

    sum_energy += e_vec[idx]*e_vec[idx];
    ++idx;
    for (int dim=0; dim<6; ++dim) {
      sum_stress += e_vec[idx]*e_vec[idx];
      ++idx;
    }
    for (int i=0; i<cell.natoms; ++i) {
      for (int dim=0; dim<3; ++dim) {
        sum_forces += e_vec[idx]*e_vec[idx];
        ++idx;
      }
    }
  }

  // Rest of punishment errors
  for (int i=idx; i<esize; ++i)
    sum_punish += e_vec[i]*e_vec[i];

  double e_sum = sum_energy + sum_stress + sum_forces + sum_punish;

  int npunish = esize - 3*natoms - 6*ncells - ncells;

  std::ofstream error_fs (oss.str(), std::ios_base::out);
  error_fs << "Total error sum " << std::fixed << e_sum << " count " << esize << " ("
           << 3*natoms << " forces, " << 6*ncells << " stresses, " << ncells << " energies";
  if (npunish == 0)
    error_fs << ")" << std::endl;
  else
    error_fs << ", " << npunish << " punishments)" << std::endl;

  error_fs << "Sum of force errors  = " << std::fixed << sum_forces << std::endl;
  error_fs << "Sum of energy errors = " << std::fixed << sum_energy << std::endl;
  error_fs << "Sum of stress errors = " << std::fixed << sum_stress << std::endl;
  if (npunish)
    error_fs << "Sum of punish errors = " << std::fixed << sum_punish << std::endl;

  error_fs.close();

  return;
}

/* ----------------------------------------------------------------------
   write punishment/constraint errors (none for this pot type)
------------------------------------------------------------------------- */

void Potential::write_punish(const ErrorVec& e_vec, String filename, int num) const
{
  return;
}

/* ----------------------------------------------------------------------
  write out extra data (none for this pot type)
------------------------------------------------------------------------- */

void Potential::write_extras(const Comm& comm, String filename, int num)
{
  return;
}

/* ----------------------------------------------------------------------
   write out actual potential
------------------------------------------------------------------------- */

std::ostream& Potential::write_pots(std::ostream& oss) const
{
  oss << "#P" << std::endl;

  for (PotFns pot_fn : pot_fns_)
    pot_fn.write_pot_fns(oss);

  oss << std::endl;

  return oss;
}

/* ----------------------------------------------------------------------
   write out potential in LAMMPS format
------------------------------------------------------------------------- */

std::ostream& Potential::write_lmps(std::ostream& oss) const
{
  for (PotFns pot_fn : pot_fns_)
    pot_fn.write_lmp_fns(oss);

  oss << std::endl;

  return oss;
}

/* ----------------------------------------------------------------------
   iostream manipulator for output flag
------------------------------------------------------------------------- */

std::ostream& MEAMZ_NS::operator<<(std::ostream& stream, const Potential::Output flag)
{
  Potential::set_output_flag(flag);
  return stream;
}

/* ----------------------------------------------------------------------
   output potential
------------------------------------------------------------------------- */

std::ostream& MEAMZ_NS::operator<<(std::ostream& stream, const Potential& pot)
{
  return pot.write(stream);
}
