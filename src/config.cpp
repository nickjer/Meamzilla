#include "config.h"
#include "input.h"
#include "universe.h"
#include "error.h"
#include "pot_list.h"
#include "style_atom.h"
#include "style_pair.h"
#include "style_triplet.h"

#include <iterator>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Config::Config(class Meamzilla *mmz) : Pointers(mmz), total_natoms(0), ncells(0), config_file_(""),
                                       energy_weight_(0.0), stress_weight_(0.0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

Config::~Config()
{
  // Atoms are deleted when the cell is deleted
  for (Cell*& cell_ptr : cells)
    if (cell_ptr) delete cell_ptr;
}

/* ----------------------------------------------------------------------
   initialize config
------------------------------------------------------------------------- */

void Config::init()
{
  if (mmz->input) {
    mmz->input->parse("config",1,config_file_);
    mmz->input->parse("energy_weight",0,energy_weight_);
    mmz->input->parse("stress_weight",0,stress_weight_);
  }

  // Read config file
  StringList conf_lines = mmz->universe->read_file(config_file_);

  // Make an atom pointer template depending on pot_type
  Atom *atom_template = nullptr;
  if (0) return;
#define ATOM_CLASS
#define AtomStyle(key,Class) \
  else if (mmz->potlist->get_pot_type() == #key) atom_template = new Class();
#include "style_atom.h"
#undef AtomStyle
#undef ATOM_CLASS
  else {
    if (mmz->universe->is_root())
      std::cerr << "ERROR: Invalid atom style - " << mmz->potlist->get_pot_type() << std::endl;
    mmz->universe->abort_all();
  }

  // Parse config file
  int line_num = 0;
  int last_line = conf_lines.size();
  do {
    try {

      // Try to read in the cell
      Cell *cell = new Cell;
      line_num = cell->read_cell(conf_lines, line_num, mmz->potlist->get_ntypes(), atom_template);
      cell->find_scale(mmz->potlist->get_max_rcut()); // find scale of cell that holds radial cutoff
      push_back(cell); // add this cell to list

    } catch(Error& ex) {

      // If it fails output the error
      if (mmz->universe->is_root())
        std::cerr << ex.what() << "cell=" << ncells << " (file=" << config_file_ << ")" << std::endl;
      mmz->universe->abort_all();

    }
  } while (line_num < last_line);

  delete atom_template; // clean up

  if (mmz->universe->is_root())
    std::cout << std::endl << "Total # atoms = " << total_natoms << " inside a total # cells = " << ncells << std::endl;

  // Find pairs for atoms
  find_pairs(atoms);

  // Find triplets for each atom in supercell if triplet is needed for potential
  if (0) return;
#define TRIPLET_CLASS
#define TripletStyle(key,Class) \
  else if (mmz->potlist->get_pot_type() == #key) find_triplets(atoms, 0);
#include "style_triplet.h"
#undef TripletStyle
#undef TRIPLET_CLASS

  // Delete all the pairs and triplets. In order to save space
  // we only need the pairs and triplets for the atoms on our proc.
  // We also only need npairs and ntriplets for sorting atoms on each proc.
  // We delete all pairs/triplets first to avoid memory fragmentation. We
  // find pairs and triplets later in atom_vec.init() call.
  delete_data(atoms);

  return;
}

/* ----------------------------------------------------------------------
   get energy weight
------------------------------------------------------------------------- */

double Config::get_energy_weight() const
{
  return energy_weight_;
}

/* ----------------------------------------------------------------------
   get stress weight
------------------------------------------------------------------------- */

double Config::get_stress_weight() const
{
  return stress_weight_;
}

/* ----------------------------------------------------------------------
   find pairs for given list of atoms
------------------------------------------------------------------------- */

void Config::find_pairs(std::vector<Atom *>& atom_list)
{
  // Make a pair pointer template depending on pot_type
  Pair *pair_template = nullptr;
  if (0) return;
#define PAIR_CLASS
#define PairStyle(key,Class) \
  else if (mmz->potlist->get_pot_type() == #key) pair_template = new Class();
#include "style_pair.h"
#undef PairStyle
#undef PAIR_CLASS
  else {
    if (mmz->universe->is_root())
      std::cerr << "ERROR: Invalid pair style - " << mmz->potlist->get_pot_type() << std::endl;
    mmz->universe->abort_all();
  }

  double max_rcut = mmz->potlist->get_max_rcut();

  for (Atom*& atom_ptr : atom_list) {   // go through list of atoms computing pairs
    Atom &atom = *atom_ptr;
    atom.npairs = 0;
    Cell &cell = *cells[atom.cell_idx];

    // Get scaling factors
    Vect cell_scale = cell.cell_scale;

    // Get lattice vectors
    std::vector<Vect> a = cell.a;

    // Check that nothing crazy happens
    if ( cell_scale.x <= 0 || cell_scale.y <= 0 || cell_scale.z <= 0 ) {
      if (mmz->universe->is_root())
        std::cerr << "ERROR: Cell #" << cell.cell_idx << " is too small or potential is too big" << std::endl;
      mmz->universe->abort_all();
    }

    // Do we double count pairs for this potential
    // set to 0 if we do, otherwise set to atom index in the cell
    int j0;
    if (mmz->potlist->dbl_cnt_pair()) j0 = 0;
    else j0 = atom.atom_idx;   // position of atom in cell

    for (int j=j0; j<cell.natoms; ++j) {
      for ( int ix = -cell_scale.x; ix <= cell_scale.x; ++ix ) {
        for ( int iy = -cell_scale.y; iy <= cell_scale.y; ++iy ) {
          for ( int iz = -cell_scale.z; iz <= cell_scale.z; ++iz ) {

            // Make sure atom_j isn't original atom_i
            if ((atom.atom_idx == j) && (ix == 0) && (iy == 0) && (iz == 0)) continue;

            Vect dist;  // Distance vector between atoms
            dist.x = cell.atoms[j]->pos.x + ix*a[0].x + iy*a[1].x + iz*a[2].x - atom.pos.x;
            dist.y = cell.atoms[j]->pos.y + ix*a[0].y + iy*a[1].y + iz*a[2].y - atom.pos.y;
            dist.z = cell.atoms[j]->pos.z + ix*a[0].z + iy*a[1].z + iz*a[2].z - atom.pos.z;

            double r = mag(dist);

            // First check that this pair lies inside maximum radius
            if (r >= max_rcut) continue;

            // Make temporary pair
            Pair *pair_ptr = pair_template->clone();
            pair_ptr->setup_pair(cell.atoms[j],dist);

            // Check that this pair is within boundaries of potentials
            // and initialize its properties
            try {
              if (pair_ptr->check_pair(atom_ptr, mmz->potlist->pot_template()))
                atom.push_back(pair_ptr);
              else
                delete pair_ptr;
            } catch(Error& ex) {
              // If it fails output the error
              if (mmz->universe->is_root()) {
                std::cerr << ex.what() << "cell=" << cell.cell_idx;
                std::cerr << " atom=" << atom.atom_idx;
                std::cerr << " neigh=" << j;
                std::cerr << " (file=" << config_file_ << ")" << std::endl;
              }
              mmz->universe->abort_all();
            }
          } // cell scale in z-direction
        } // cell scale in y-direction
      } // cell scale in x-direction
    } // loop over neighbor atoms
  } // loop over atom_list

  delete pair_template;   // clean up

  return;
}

/* ----------------------------------------------------------------------
   find triplets for atoms and whether to save data
------------------------------------------------------------------------- */

void Config::find_triplets(std::vector<Atom *>& atom_list, int save_data)
{
  // Make a triplet pointer template depending on pot_type
  Triplet *triplet_template = nullptr;
  if (0) return;
#define TRIPLET_CLASS
#define TripletStyle(key,Class) \
  else if (mmz->potlist->get_pot_type() == #key) triplet_template = new Class();
#include "style_triplet.h"
#undef TripletStyle
#undef TRIPLET_CLASS
  else {
    if (mmz->universe->is_root())
      std::cerr << "ERROR: Invalid triplet style - " << mmz->potlist->get_pot_type() << std::endl;
    mmz->universe->abort_all();
  }

  for (Atom*& atom_ptr : atom_list) {   // go through list of atoms computing pairs
    Atom &atom = *atom_ptr;
    atom.ntriplets = 0;

    // Now loop through all combinations of neighbors for this atom
    for (int j=0; j<atom.npairs-1; ++j) {
      for (int k=j+1; k<atom.npairs; ++k) {

        // Make temporary triplet
        Triplet *triplet_ptr = triplet_template->clone();
        triplet_ptr->setup_triplet(atom.pairs[j], atom.pairs[k]);

        // Check that this is a triplet and initialize its properties
        if (triplet_ptr->check_triplet(atom_ptr, mmz->potlist->pot_template())) {
          if (save_data) {
            atom.push_back(triplet_ptr);
          } else {
            delete triplet_ptr;
            ++atom.ntriplets; // count this triplet but don't save it
          }
        }
        else
          delete triplet_ptr;
      }
    }
  }

  return;
}

/* ----------------------------------------------------------------------
  delete memory used by pairs and triplets
------------------------------------------------------------------------- */

void Config::delete_data(std::vector<Atom *>& atom_list)
{
  for (Atom*& atom_ptr : atom_list) {   // go through list of atoms deleting pairs/triplets
    for (Triplet*& triplet_ptr : atom_ptr->triplets) delete triplet_ptr;
    atom_ptr->triplets.clear();

    for (Pair*& pair_ptr : atom_ptr->pairs) delete pair_ptr;
    atom_ptr->pairs.clear();
  }

  return;
}

/* ----------------------------------------------------------------------
   add cell to list of cells
------------------------------------------------------------------------- */

void Config::push_back(Cell *cell)
{
  // update cell with cell idx and # of atoms up to this point
  cell->update_cell(ncells, total_natoms);

  cells.push_back(cell);         // add cell to list
  ++ncells;                      // count total number of cells
  total_natoms += cell->natoms;   // count total number of atoms
  // Add atoms from cell to list of atoms
  for (Atom*& atom_ptr : cell->atoms)
    atoms.push_back(atom_ptr);

  return;
}

/* ----------------------------------------------------------------------
   remove specified cell from list
------------------------------------------------------------------------- */

void Config::erase_cell(int idx)
{
  if (cells[idx]) delete cells[idx];  // this deletes the cell and all of its atoms
  cells.erase(cells.begin() + idx);   // remove this cell from the list of cells
  --ncells;  // keep count of cells

  // Erase empty atom pointers from atoms list and keep track of total_natoms
  std::vector<Atom *>::iterator erase_begin = atoms.end(), erase_end = atoms.end();
  for (std::vector<Atom *>::iterator aiter = atoms.begin(); aiter != atoms.end(); ++aiter) {

    if (!*aiter && erase_begin == atoms.end()) erase_begin = aiter; // set position of erasure beginning

    if (*aiter && erase_begin != atoms.end()) erase_end = aiter;    // set position of erasure end
  }
  atoms.erase(erase_begin, erase_end);               // erase these atoms from list of atoms

  // Need to update cell_idx and global_atom_idx for cells and atoms after the deletion point
  // Don't forget to get updated number of atoms
  total_natoms = 0;
  for (int i=0; i<idx; ++i) total_natoms += cells[i]->natoms;  // get # of atoms up to deletion point
  for (int i=idx; i<ncells; ++i) {
    cells[i]->update_cell(i, total_natoms); // update cell with cell idx and # of atoms up to this point
    total_natoms += cells[i]->natoms;    // count total number of atoms
  }

  return;
}
