/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#ifndef MEAMZ_CELL_H
#define MEAMZ_CELL_H

#include "vect.h"
#include "stress.h"
#include "atom.h"
#include <string>
#include <vector>

namespace MEAMZ_NS
{

typedef std::string String;
typedef std::vector<String> StringList;

class Cell
{
public:
  int cell_idx;                  // Index of supercell in config

  int natoms;                    // Number of atoms in cell
  std::vector<Atom *> atoms;     // List of atoms that make up cell
  std::vector<Vect> a;           // 3-D lattice vectors of cell

  double vol;                    // Volume of cell
  double weight;                 // Weight of cell in fit
  double eweight;                // Weight of cell energy in fit
  double sweight;                // Weight of cell stress in fit
  double energy;                 // Energy of cell
  double energy0;                // User supplied energy of cell
  Stress stress;                 // Stress tensor
  Stress stress0;                // User supplied stress tensor

  int useforce;                  // Whether we use the forces supplied in the fitting

  Vect cell_scale;               // How large to scale the cell in 3 dimensions to hold radial potential

  Cell();
  virtual ~Cell();

  int read_cell(const StringList&, int, int, Atom *);       // Read in a cell from a configuration file
                                                            //    Could be made virtual to read in different formats
  void update_cell(int, int);         // Update cell and atoms with cell_idx and global_atom_idx
  void find_scale(double);            // Compute the scale of the cell that holds the radial potential

protected:

private:
  void compute_vol();           // Compute volume of supercell

};

}; // namespace MEAMZ_NS

#endif // MEAMZ_CELL_H
