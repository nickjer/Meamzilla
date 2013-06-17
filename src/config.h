#ifndef MEAMZ_CONFIG_H
#define MEAMZ_CONFIG_H

#include "pointers.h"
#include "cell.h"
#include "atom.h"

namespace MEAMZ_NS
{

class Config : protected Pointers
{
public:
  int total_natoms;             // Number of atoms in the configuration database
  std::vector<Atom *> atoms;    // List of atoms

  int ncells;                   // Number of cells in the configuration database
  std::vector<Cell *> cells;    // List of cells

  Config(class Meamzilla *);
  virtual ~Config();

  void init();                      // Initialize config

  double get_energy_weight() const; // Get energy weight
  double get_stress_weight() const; // Get stress weight

  void find_pairs(std::vector<Atom *>&);              // Find pairs for given list of atoms
  void find_triplets(std::vector<Atom *>&, int = 0);  // Find triplets for atoms and whether to save data
  void delete_data(std::vector<Atom *>&);             // Delete memory used by pairs and triplets

protected:

private:
  String config_file_;    // Name of config file
  double energy_weight_;  // Energy weight of whole configuration
  double stress_weight_;  // Stress weight of whole configuration

  void push_back(Cell *);       // Add cell to the end of the list of cells
  void erase_cell(int);         // Remove specified cell from list
};

} // namespace MEAMZ_NS

#endif // MEAMZ_CONFIG_H
