/* ----------------------------------------------------------------------
   MEAMZ - MEAM optimiZer
   http://github.com/nickjer/Meamzilla
   Jeremy Nicklas, nicklas.2@buckeyemail.osu.edu

   Copyright (2013) Jeremy Nicklas.  All rights reserved.  This code
   is not free to use without the express permission of the author.

   See the README file in the top-level MEAMZ directory.
------------------------------------------------------------------------- */

#include "atom_vec.h"
#include "config.h"
#include "universe.h"
#include "pot_list.h"

#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

AtomVec::AtomVec(class Meamzilla *mmz) : Pointers(mmz), natoms(0), npairs_(0), ntriplets_(0), nclusters_(0)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

AtomVec::~AtomVec()
{
  //dtor
}

/* ----------------------------------------------------------------------
   initialize atom vector
------------------------------------------------------------------------- */

void AtomVec::init()
{
  // First create vector of all atoms
  AtomVec tmp_all_atoms(mmz);
  for (Atom*& atom_ptr : mmz->config->atoms)
    tmp_all_atoms.push_back(atom_ptr);

  // Now we need to sort this vector in descending order based on
  // user's choice of atom's potential (for example: # pairs or # triplets)
  std::sort(tmp_all_atoms.atoms.begin(),tmp_all_atoms.atoms.end(),
            [](const Atom *la, const Atom *ra){
              return la->cmp_value() > ra->cmp_value();
            });

  // Ideal MPI system would have this many clusters per processor
  int nprocs = mmz->universe->comm_inter.get_size();
  double clust_per_proc = tmp_all_atoms.nclusters_*1.0/nprocs;

  // Make sure there isn't too many procs
  if (nprocs>tmp_all_atoms.natoms) {
    if (mmz->universe->comm_inter.is_root())
      std::cerr << "ERROR: More processors than atoms" << std::endl;
    mmz->universe->comm_inter.abort_all();
  }

  // Make a clean list of atoms and its number of clusters for each proc
  std::vector<AtomVec> proc_atomvecs(nprocs, AtomVec(mmz));

  for (AtomVec& proc_atomvec : proc_atomvecs) {
    // Store first atom in list for this proc
    proc_atomvec.push_back(tmp_all_atoms.atoms[0]);
    std::vector<int> used_idx(1,0); // keep track of atoms already used for future deletion

    // Cycle through each atom adding it to list of atoms for proc if
    // the sum of clusters doesn't surpass the ideal # of clusters per proc
    for (int j=1; j<tmp_all_atoms.natoms; ++j) {
      if (proc_atomvec.nclusters_+tmp_all_atoms.atoms[j]->cmp_value()>clust_per_proc)
        continue;
      proc_atomvec.push_back(tmp_all_atoms.atoms[j]);
      used_idx.push_back(j);
    }

    // Erase elements we used
    for (int j=0; j<int(used_idx.size()); ++j) tmp_all_atoms.erase(used_idx[j]-j);
    used_idx.clear();
  }

  // Sort the list of atoms for each proc based on # of clusters: smallest to largest
  std::sort(proc_atomvecs.begin(),proc_atomvecs.end());

  // Any left over atoms we just append them to each proc list (smaller ones first)
  for (int i=0; i<tmp_all_atoms.natoms; ++i)
    proc_atomvecs[i%nprocs].push_back(tmp_all_atoms.atoms[i]);

  // Now assign an atom list to each processor
  int me = mmz->universe->comm_inter.get_rank();
  *this = proc_atomvecs[me];

  if (mmz->universe->comm_intra) {
    for (int group=0; group<mmz->universe->get_ngroups(); ++group) {
      if (group == mmz->universe->get_group()) {
        std::cout << std::endl;
        std::cout << "Group #" << group << ":" << std::endl;
        std::cout << std::setw(6)  << std::right << "Node";
        std::cout << std::setw(12) << std::right << "# Atoms";
        std::cout << std::setw(12) << std::right << "# Pairs";
        std::cout << std::setw(15) << std::right << "# Triplets" << std::endl;
        std::cout << std::setfill('-') << std::setw(46) << std::left << " ";
        std::cout << std::setfill(' ');
        std::cout << std::endl;

        for (int i=0; i<nprocs; ++i) {
          std::cout << std::setw(6)  << std::right << i;
          std::cout << std::setw(12) << std::right << proc_atomvecs[i].natoms;
          std::cout << std::setw(12) << std::right << proc_atomvecs[i].npairs_;
          std::cout << std::setw(15) << std::right << proc_atomvecs[i].ntriplets_;
          std::cout << std::endl;
        }

        std::cout << std::endl;
      }
      mmz->universe->comm_intra.barrier();
    }
  }
  mmz->universe->comm_all.barrier();

  // Find pairs and triplets again
  mmz->config->find_pairs(atoms);

  // Find triplets for each atom in supercell if triplet is needed for potential
  if (0) return;
#define TRIPLET_CLASS
#define TripletStyle(key,Class) \
  else if (mmz->potlist->get_pot_type() == #key) mmz->config->find_triplets(atoms, 1);
#include "style_triplet.h"
#undef TripletStyle
#undef TRIPLET_CLASS

  return;
}

/* ----------------------------------------------------------------------
   adds atom to list of atoms in this object
------------------------------------------------------------------------- */

void AtomVec::push_back(Atom *atom)
{
  atoms.push_back(atom);
  ++natoms;
  npairs_ += atom->npairs;
  ntriplets_ += atom->ntriplets;
  nclusters_ += atom->cmp_value();
}

/* ----------------------------------------------------------------------
   erases atom from list of atoms in this object
------------------------------------------------------------------------- */

void AtomVec::erase(int idx)
{
  --natoms;
  npairs_ -= atoms[idx]->npairs;
  ntriplets_ -= atoms[idx]->ntriplets;
  nclusters_ -= atoms[idx]->cmp_value();
  atoms.erase(atoms.begin() + idx);
}
