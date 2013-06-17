#include "opt_list.h"
#include "input.h"
#include "universe.h"
#include "pot_list.h"
#include "style_opt.h"

#include <sstream>

using namespace MEAMZ_NS;

/* ---------------------------------------------------------------------- */

OptList::OptList(class Meamzilla *mmz) : Pointers(mmz)
{
  //ctor
}

/* ---------------------------------------------------------------------- */

OptList::~OptList()
{
  for (Optimizer*& opt : opts_) if (opt) delete opt;
}

/* ----------------------------------------------------------------------
   initialize list of optimizers
------------------------------------------------------------------------- */

void OptList::init()
{
  // Get optstyle line from input file (default: none)
  StringList optstyle_list;
  if (mmz->input) {
    mmz->input->parse_stringlist("optstyle", 0, optstyle_list);
  }

  // Create optimizer class
  nopts_ = optstyle_list.size();
  for (int i=0; i<nopts_; ++i) {
    if (0) return;
#define OPTIMIZER_CLASS
#define OptimizerStyle(key,Class) \
    else if (optstyle_list[i] == #key) opts_.push_back(new Class(mmz, &std::cout));
#include "style_opt.h"
#undef OptimizerStyle
#undef OPTIMIZER_CLASS
    else {
      if (mmz->universe->is_root())
        std::cerr << "ERROR: Invalid optimizer style - " << optstyle_list[i] << std::endl;
      mmz->universe->abort_all();
    }

    opts_.back()->init();
  }

  return;
}

/* ----------------------------------------------------------------------
   compute list of optimizers using list of potentials
------------------------------------------------------------------------- */

void OptList::compute(PotList& pot_list)
{
  for (Optimizer*& opt : opts_)
    opt->compute(pot_list);

  return;
}
