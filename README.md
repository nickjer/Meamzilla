Meamzilla
=========

Molecular dynamics potential fitting code. Fits a classical potential to a configuration database of ab-initio forces, energies, and stresses.

Requirements
============

gcc 4.7.2+ (due to C++11 standards in the code)
MPI
MKL/ACML/LibSci

Installation
============

Edit the Makefile with proper C++ compiler and libraries. There are a few example Makefile's provided for a variety of supercomputers. Then type:

make -j <# cores>
