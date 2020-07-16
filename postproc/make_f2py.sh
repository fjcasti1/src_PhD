#!/bin/bash
# ======================================================================
# Wraps fortran 90 code for python
#
# Environment Requires:
#   module load intel/2019.4
#   and appropriate python env
# ======================================================================
#MODULE="intel/2019.4"
#
#if lsmod | grep "$MODULE" &> /dev/null ; then
#  echo "Module $MODULE is loaded"
#else
#  echo "Loading module $MODULE"
##  module load intel/2019.4
#fi

module purge
module load intel/2019.4

prefix=reader
prefix=myParaview2
opts=(
  ${prefix}.F90     # Fortran file to compile
  -m ${prefix}      # Resultant module file
  -c                # Allow f2py to accept compiler options
  --fcompiler=intel # Specify compilier family (intel)
  # Specify optimization level, allow vectorization and parallel mkl
  # NOTE: MKL_NUM_THREADS env variable will limit parallel execution
  --opt="-O3 -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512 -mkl=parallel"
)
f2py "${opts[@]}"

# To load paloma's
# ifort -03 -mkl=parallel -xAVX2 pv_paloma.f90 -o pv_paloma
#
# To run it
# ./pv_paloma < pv_paloma_input
