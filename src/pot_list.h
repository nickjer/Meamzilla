#ifndef MEAMZ_POT_LIST_H
#define MEAMZ_POT_LIST_H

#include "pointers.h"
#include "potential.h"

namespace MEAMZ_NS
{

class PotList : protected Pointers
{
public:
  PotList(class Meamzilla *);
  PotList(const PotList&);
  virtual ~PotList();

  PotList& operator=(const PotList&);

  void init();                  // Initialize this class
  int get_ntypes() const;       // Get number of atom types (useful for alloys)
  double get_max_rcut() const;  // Get maximum radial cutoff from all potentials
  int get_rescale() const;      // Get rescale value
  String get_pot_type() const;  // Get potential type for list of potentials
  String get_tmp_file() const;  // Get name of temporary pot output file
  String get_end_file() const;  // Get name of ending pot file for output
  String get_dat_file() const;  // Get name of data file
  String get_lmp_file() const;  // Get name of lammps output file
  int dbl_cnt_pair() const;     // Decide whether we should double count pairs based on potential type chosen

  int get_npots() const;        // Get total number of potentials in list
  int get_group_npots() const;  // Get number of potentials this group sees
  int get_group(int) const;     // Get group number corresponding to global potential idx

  void write_pots(String);      // Output list of pots to a user specified file
  void write_lmps(String);      // Output pots in lammps format to user specified file
  void write_data(String);      // Output physics data to a user specified file
  void sort();                  // Sort list of potentials based on error sum

  Potential * pot_template() const;   // Get first potential as template for config pairs/triplets
  Potential *& at(int);               // Give pointer to potential depending on group processor is in
  Potential *& operator[](int);       // Give pointer to potential depending on group processor is in
  Potential *& global_at(int);        // Give pointer to potential using global idx
  void resize(int);                   // Resize the list of pots
  void push_back(Potential *);        // Add pot to list of pots

protected:

private:
  int ntypes_;                      // Number of atom types used for alloys
  int dbl_cnt_pair_;                // Integer specifying if potentials rely on double counting of pairs or not
  int rescale_;                     // Integer specifying if potentials should be rescaled

  String pot_file_;                 // Starting pot file
  String tmp_file_;                 // Temporary pot output file
  String end_file_;                 // Ending pot file for output
  String dat_file_;                 // Data file
  String lmp_file_;                 // Lammps file for output

  String pot_type_;                 // Type of potential used

  String comment_line_;             // Comment line used in potential file

  int npots_;                       // Number of potentials
  std::vector<Potential *> pots_;   // List of potentials (must be same type of potential)
};

} // namespace MEAMZ_NS

#endif // MEAMZ_POT_LIST_H
