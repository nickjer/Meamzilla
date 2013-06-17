#include "cell.h"
#include "error.h"

#include <sstream>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

Cell::Cell() : cell_idx(0), natoms(0), a(3), weight(1), eweight(1), sweight(1),
               energy(0), energy0(0), useforce(1)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

Cell::~Cell()
{
  for (Atom*& atom_ptr : atoms)
    if (atom_ptr) delete atom_ptr;
}

/* ----------------------------------------------------------------------
   initialize config
------------------------------------------------------------------------- */

int Cell::read_cell(const StringList& lines, int line_num, int ntypes, Atom *atom_template)
{
  std::istringstream iss (lines[line_num], std::istringstream::in);

  int end_line = lines.size(); // last line in file

  // First line of cell must be (#N <num atoms> <whether to use forces?>)
  std::string label;
  iss >> label >> natoms >> useforce;
  if (label != "#N")
    throw Error("ERROR: number line (#N ...) is missing - ");
  if (iss.fail())
    throw Error("ERROR: number line (#N ...) is incorrect - ");
  if (line_num+1 == end_line)
    throw Error("ERROR: unexpected end of config file - ");

  // Read in header of cell, stopping at #F
  while(lines[++line_num].substr(0,2) != "#F") {

    if (line_num+1 == end_line)
      throw Error("ERROR: unexpected end of config file - ");

    iss.clear();
    iss.str(lines[line_num]);

    iss >> label;
    if (label=="#X") {
      iss >> a[0];
      if (iss.fail()) throw Error("ERROR: box_x vector (#X ...) incorrect - ");
    } else if (label=="#Y") {
      iss >> a[1];
      if (iss.fail()) throw Error("ERROR: box_y vector (#Y ...) incorrect - ");
    } else if (label=="#Z") {
      iss >> a[2];
      if (iss.fail()) throw Error("ERROR: box_z vector (#Z ...) incorrect - ");
    } else if (label=="#E") {
      iss >> energy0;
      if (iss.fail()) throw Error("ERROR: energy (#E ...) incorrect - ");
    } else if (label=="#S") {
      iss >> stress0;
      if (iss.fail()) throw Error("ERROR: stress tensor (#S ...) incorrect - ");
    } else if (label=="#W") {
      iss >> weight;
      if (iss.fail()) throw Error("ERROR: configuration weight (#W ...) incorrect - ");
    } else if (label=="#EW") {
      iss >> eweight;
      if (iss.fail()) throw Error("ERROR: configuration energy weight (#EW ...) incorrect - ");
    } else if (label=="#SW") {
      iss >> sweight;
      if (iss.fail()) throw Error("ERROR: configuration stress weight (#SW ...) incorrect - ");
    }
  }

  if (line_num+natoms >= end_line)
    throw Error("ERROR: unexpected end of config file - ");

  compute_vol();  // compute cell's volume

  atoms.resize(natoms, nullptr);  // resize list of atoms
  for (Atom*& atom_ptr : atoms)
    atom_ptr = atom_template->clone();

  // Add atoms
  for (int i=0; i<natoms; ++i) {
    String atom_line = lines[++line_num];

    if (atom_line[0] == '#')
      throw Error("ERROR: too few atoms - ");

    iss.clear();
    iss.str(atom_line);

    int typ;
    Vect pos, force;
    iss >> typ >> pos >> force;
    if (iss.fail()) {
      std::ostringstream oss;
      oss << "ERROR: atom line incorrect - atom=" << i << " ";
      throw Error(oss.str());
    }
    if (typ < 0 || typ >= ntypes) {
      std::ostringstream oss;
      oss << "ERROR: incorrect atom type specified - atom=" << i << " ";
      throw Error(oss.str());
    }

    atoms[i]->setup_atom(i,typ,pos,force);
  }

  // Peek ahead a line and make sure it starts with #N for new cell
  ++line_num;
  if (line_num < end_line)
    if (lines[line_num].substr(0,2)!="#N")
      throw Error("ERROR: too many atoms - ");
  --line_num; // after peeking ahead, go back a line for consistency

  return ++line_num;
}

/* ----------------------------------------------------------------------
   update cell and atoms with cell_idx and global_atom_idx
------------------------------------------------------------------------- */

void Cell::update_cell(int _cell_idx, int _global_atom_idx)
{
  cell_idx = _cell_idx;

  for (Atom*& atom_ptr : atoms) {
    atom_ptr->cell_idx = _cell_idx;
    atom_ptr->global_idx = _global_atom_idx++;
  }

  return;
}

/* ----------------------------------------------------------------------
   compute the scale of the cell that holds the radial potential
------------------------------------------------------------------------- */

void Cell::find_scale(double max_rcut)
{
  // Compute inverse heights to get scaling factors
  Vect inv_height;
  inv_height.x = mag(a[1] ^ a[2]) / vol;    // = |a1 x a2|/V
  inv_height.y = mag(a[2] ^ a[0]) / vol;    // = |a2 x a0|/V
  inv_height.z = mag(a[0] ^ a[1]) / vol;    // = |a0 x a1|/V

  // Compute scaling factors
  cell_scale.x = int(std::ceil(max_rcut * inv_height.x));
  cell_scale.y = int(std::ceil(max_rcut * inv_height.y));
  cell_scale.z = int(std::ceil(max_rcut * inv_height.z));
}

/* ----------------------------------------------------------------------
   compute volume of supercell
------------------------------------------------------------------------- */

void Cell::compute_vol()
{
  vol = std::abs( a[0] * ( a[1] ^ a[2] ) );   // V = |a1 . ( a2 x a3 )|
}

