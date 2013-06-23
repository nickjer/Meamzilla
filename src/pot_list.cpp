/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "pot_list.h"
#include "input.h"
#include "universe.h"
#include "error.h"
#include "config.h"
#include "style_pot.h"

#include <fstream>
#include <algorithm>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

PotList::PotList(class Meamzilla *mmz) : Pointers(mmz), ntypes_(1), dbl_cnt_pair_(0), rescale_(2),
                                         pot_file_("start"), tmp_file_("temp"), end_file_("end"),
                                         dat_file_("data"), lmp_file_(""), comment_line_("#"), npots_(0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

PotList::PotList(const PotList& pot_list) : Pointers(pot_list.mmz), ntypes_(pot_list.ntypes_),
                                            dbl_cnt_pair_(pot_list.dbl_cnt_pair_), rescale_(pot_list.rescale_),
                                            pot_file_(pot_list.pot_file_), tmp_file_(pot_list.tmp_file_),
                                            end_file_(pot_list.end_file_), dat_file_(pot_list.dat_file_),
                                            lmp_file_(pot_list.lmp_file_), pot_type_(pot_list.pot_type_),
                                            comment_line_(pot_list.comment_line_), npots_(pot_list.npots_)
{
  pots_.resize(npots_, nullptr);
  for (int i=0; i<npots_; ++i)
    pots_[i] = pot_list.pots_[i]->clone();
}

/* ---------------------------------------------------------------------- */

PotList::~PotList()
{
  for (Potential*& pot : pots_) if (pot) delete pot;
}

/* ---------------------------------------------------------------------- */

PotList& PotList::operator=(const PotList& pot_list)
{
  mmz = pot_list.mmz;
  ntypes_ = pot_list.ntypes_;
  dbl_cnt_pair_ = pot_list.dbl_cnt_pair_;
  rescale_ = pot_list.rescale_;
  pot_file_ = pot_list.pot_file_;
  tmp_file_ = pot_list.tmp_file_;
  end_file_ = pot_list.end_file_;
  dat_file_ = pot_list.dat_file_;
  lmp_file_ = pot_list.lmp_file_;
  pot_type_ = pot_list.pot_type_;
  comment_line_ = pot_list.comment_line_;
  npots_ = pot_list.npots_;

  pots_.resize(npots_, nullptr);
  for (int i=0; i<npots_; ++i) {
    if (pots_[i]) delete pots_[i];
    pots_[i] = pot_list.pots_[i]->clone();
  }

  return *this;
}

/* ----------------------------------------------------------------------
   initialize list of potentials
------------------------------------------------------------------------- */

void PotList::init()
{
  if (mmz->input) {
    mmz->input->parse("rescale",    0, rescale_);
    mmz->input->parse("startpot",   1, pot_file_);
    mmz->input->parse("tempfile",   1, tmp_file_);
    mmz->input->parse("endpot",     1, end_file_);
    mmz->input->parse("datafile",   0, dat_file_);
    mmz->input->parse("lammpsfile", 0, lmp_file_);
  }

  // Read pot file
  StringList pot_lines = mmz->universe->read_file(pot_file_);

  // Parse potential list file
  int line_num = 0;
  int last_line = pot_lines.size();
  try {
    // Read in comment line first
    comment_line_ = pot_lines[line_num++];

    // Try to read in the potential type followed by the potential
    if (pot_lines[line_num].substr(0,2)!="#T") throw Error("ERROR: potential type line (#T ...) missing - ");

    std::istringstream iss (pot_lines[line_num++], std::istringstream::in);
    std::string label;
    iss >> label >> pot_type_ >> ntypes_;
    if (iss.fail()) throw Error("ERROR: pot type line (#T ...) is incorrect - ");

    // Create potential class and check user defined types
    Potential *pot = nullptr;
    if (0) return;
#define POTENTIAL_CLASS
#define PotentialStyle(key,Class) \
    else if (pot_type_ == #key) { pot = new Class(mmz, ntypes_); pot->init(); }
#include "style_pot.h"
#undef PotentialStyle
#undef POTENTIAL_CLASS
    else {
      if (mmz->universe->is_root())
        std::cerr << "ERROR: invalid potential style \"" << pot_type_ << "\" - " << std::endl;
      mmz->universe->abort_all();
    }

    // Decide if this potential needs double counting of pairs or not and store that info for later use
    dbl_cnt_pair_ = pot->dbl_cnt_pair();

    // Cycle through reading each potential from list
    do {
      if (pot_lines[line_num++].substr(0,2)!="#P") throw Error("ERROR: potential line (#P ...) missing - ");

      // Create potential class
      Potential *tmp_pot = pot->clone();

      // Read in the potential
      line_num = tmp_pot->read_pot(pot_lines, line_num);

      push_back(tmp_pot); // add this pot to list
    } while (line_num < last_line);

    delete pot; // clean up

  } catch(Error& ex) {

    // If it fails output the error
    if (mmz->universe->is_root())
      std::cerr << ex.what() << "pot=" << npots_ << " (file=" << pot_file_ << ")" << std::endl;
    mmz->universe->abort_all();

  }

  // Don't allow more groups than potentials
  int ngroups = mmz->universe->get_ngroups();
  if (ngroups > npots_) {
    if (mmz->universe->is_root())
      std::cerr << "ERROR: User defined more groups (ngroups=" << ngroups << ") than potentials (npots=" << npots_ << ")" << std::endl;
    mmz->universe->abort_all();
  }

  // Check if potentials are similar enough to run calculations on them
  // using the specified pairs/triplets
  for (int i=1; i<npots_; ++i) {
    try {
      pots_[0]->check_if_similar(pots_[i]); // error will be thrown if not similar
    } catch(Error& ex) {
      // If it fails output the error
      if (mmz->universe->is_root())
        std::cerr << ex.what() << "pota=0 potb=" << i << " (file=" << pot_file_ << ")" << std::endl;
      mmz->universe->abort_all();
    }
  }

  return;
}

/* ----------------------------------------------------------------------
   get number of atom types (useful for alloys)
------------------------------------------------------------------------- */

int PotList::get_ntypes() const
{
  return ntypes_;
}

/* ----------------------------------------------------------------------
   get maximum radial cutoff from all potentials
------------------------------------------------------------------------- */

double PotList::get_max_rcut() const
{
  double max_rcut = 0.0;

  for (Potential *pot : pots_)
    max_rcut = std::max(max_rcut, pot->get_max_rcut());

  return max_rcut;
}

/* ----------------------------------------------------------------------
   get rescale value
      0 - never rescale
      1 - don't rescale in beginning / rescale throughout
      2 - rescale if necessary in beginning / rescale throughout
      3 - must rescale beginning / rescale throughout
------------------------------------------------------------------------- */

int PotList::get_rescale() const
{
  return rescale_;
}

/* ----------------------------------------------------------------------
   get potential type for list of potentials
------------------------------------------------------------------------- */

String PotList::get_pot_type() const
{
  return pot_type_;
}

/* ----------------------------------------------------------------------
   get name of temporary pot output file
------------------------------------------------------------------------- */

String PotList::get_tmp_file() const
{
  return tmp_file_;
}

/* ----------------------------------------------------------------------
   get name of ending pot file for output
------------------------------------------------------------------------- */

String PotList::get_end_file() const
{
  return end_file_;
}

/* ----------------------------------------------------------------------
   get name of data file
------------------------------------------------------------------------- */

String PotList::get_dat_file() const
{
  return dat_file_;
}

/* ----------------------------------------------------------------------
   get name of lammps output file
------------------------------------------------------------------------- */

String PotList::get_lmp_file() const
{
  return lmp_file_;
}

/* ----------------------------------------------------------------------
   decide whether we should double count pairs based on
   potential type chosen
------------------------------------------------------------------------- */

int PotList::dbl_cnt_pair() const
{
  return dbl_cnt_pair_;
}

/* ----------------------------------------------------------------------
   get total number of potentials in list
------------------------------------------------------------------------- */

int PotList::get_npots() const
{
  return npots_;
}

/* ----------------------------------------------------------------------
   get number of potentials this group sees
------------------------------------------------------------------------- */

int PotList::get_group_npots() const
{
  int group = mmz->universe->get_group();
  int ngroups = mmz->universe->get_ngroups();
  int group_npots = npots_ / ngroups;
  int remainder = npots_ % ngroups;
  if (group < remainder) ++group_npots;
  return group_npots;
}

/* ----------------------------------------------------------------------
   get group number corresponding to global potential idx
------------------------------------------------------------------------- */

int PotList::get_group(int global_idx) const
{
  int ngroups = mmz->universe->get_ngroups();
  return global_idx % ngroups;
}

/* ----------------------------------------------------------------------
   output list of pots to a user specified file
------------------------------------------------------------------------- */

void PotList::write_pots(String filename)
{
  // Communicate potentials from respective groups to all root procs (comm_intra)
  for (int i=0; i<npots_; ++i)
    pots_[i]->communicate(mmz->universe->comm_intra, get_group(i));

  if (mmz->universe->is_root()) {
    std::ofstream ofs (filename, std::ios_base::out);

    if (ofs) {
      ofs << comment_line_ << std::endl;
      ofs << "#T " << pot_type_ << " " << ntypes_ << std::endl << std::endl;
      for (Potential *pot : pots_)
        ofs << Potential::Output::Pots << *pot;
    }

    ofs.close();
  }

  return;
}

/* ----------------------------------------------------------------------
   output pots in lammps format to user specified file
------------------------------------------------------------------------- */

void PotList::write_lmps(String filename)
{
  // Communicate potentials from respective groups to all root procs (comm_intra)
  for (int i=0; i<npots_; ++i)
    pots_[i]->communicate(mmz->universe->comm_intra, get_group(i));

  if (mmz->universe->is_root()) {
    for (int i=0; i<npots_; ++i) {
      std::ostringstream oss;
      oss << filename;
      if (npots_ > 1)
      oss << "." << i;

      std::ofstream ofs (oss.str(), std::ios_base::out);

      if (ofs) {
        ofs << comment_line_ << std::endl;
        ofs << pot_type_ << std::endl;
        ofs << Potential::Output::Lmps << *pots_[i];
      }

      ofs.close();
    }
  }

  mmz->universe->comm_all.barrier();

  return;
}

/* ----------------------------------------------------------------------
   output physics data to a user specified file
------------------------------------------------------------------------- */

void PotList::write_data(String filename)
{
  if (npots_ == 1)
    pots_[0]->write_data(mmz->universe->comm_inter, filename);
  else
    for (Potential *pot : pots_)
      pot->write_data(mmz->universe->comm_inter, filename, pot->get_global_idx());

  mmz->universe->comm_all.barrier();

  return;
}

/* ----------------------------------------------------------------------
   sort list of potentials based on error sum
------------------------------------------------------------------------- */

void PotList::sort()
{
  // Communicate potentials from respective groups to all group root procs (comm_intra)
  for (int i=0; i<npots_; ++i)
    pots_[i]->communicate(mmz->universe->comm_intra, get_group(i));

  // Communicate potentials from group root proc to all procs in group
  //for (int i=0; i<npots_; ++i)
  //  pots_[i]->communicate(mmz->universe->comm_inter, mmz->universe->comm_inter.get_root());

  std::sort(pots_.begin(), pots_.end(), [](Potential *lhs, Potential *rhs){return *lhs < *rhs;});
  return;
}

/* ----------------------------------------------------------------------
   get first potential as template for config pairs/triplets
------------------------------------------------------------------------- */

Potential * PotList::pot_template() const
{
  return pots_[0];
}

/* ----------------------------------------------------------------------
   give pointer to potential depending on group processor is in
   Ex: 6 groups - 14 potentials
        Group #     Global_Idx    Local_Idx    Group_Npots
        ---------------------------------------------------
           0         0,  1,  2     0, 1, 2          3
           1         3,  4,  5     0, 1, 2          3
           2         6,  7         0, 1             2
           3         8,  9         0, 1             2
           4        10, 11         0, 1             2
           5        12, 13         0, 1             2
------------------------------------------------------------------------- */

Potential *& PotList::at(int local_idx)
{
  if (local_idx >= get_group_npots()) {
    std::cerr << "ERROR: Went outside the bounds of the potential list" << std::endl;
    mmz->universe->comm_all.abort_one();
  }
  int group = mmz->universe->get_group();
  int ngroups = mmz->universe->get_ngroups();
  int global_idx = local_idx * ngroups + group;
  return pots_[global_idx];
}

/* ----------------------------------------------------------------------
   same as at(local_idx) function except using brackets
------------------------------------------------------------------------- */

Potential *& PotList::operator[](int local_idx)
{
  return at(local_idx);
}

/* ----------------------------------------------------------------------
   give pointer to potential using global idx
------------------------------------------------------------------------- */

Potential *& PotList::global_at(int global_idx)
{
  return pots_[global_idx];
}

/* ----------------------------------------------------------------------
   resize the list of pots
------------------------------------------------------------------------- */

void PotList::resize(int size)
{
  if (npots_ < size) {
    for (int i=0; i<(size-npots_); ++i) {
      Potential *pot = pot_template()->clone();     // clone first potential in pot list
      pot->set_global_idx(npots_);                  // record position in pot list in this potential
      pots_.push_back(pot);                         // add pot to pot list
      ++npots_;
    }
  } else if (npots_ > size) {
    for (int i=size; i<npots_; ++i)
      if (pots_[i]) delete pots_[i];
    npots_ = size;
    pots_.resize(npots_);
  }

  return;
}

/* ----------------------------------------------------------------------
   add pot to list of pots
------------------------------------------------------------------------- */

void PotList::push_back(Potential *pot)
{
  pot->set_global_idx(npots_);    // record position in pot list in this potential
  pots_.push_back(pot);           // add pot to list
  ++npots_;                       // count total number of pots

  return;
}
