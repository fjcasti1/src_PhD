#!/usr/bin/env bash
job_rec="${1:?'JOB RECORD NOT SUPPLIED'}"
sourceFile="${2:?'SOURCE FILE NOT SUPPLIED'}"
timeLimit="${3:?'TIME LIMIT NOT SUPPLIED'}"
partition="${4:?'PARTITION NOT SUPPLIED'}"
MODE="${5:?'MODE NOT SUPPLIED'}"
job_com="${6:-main}"

if [ $partition == "debug" ]; then
  echo "Partition = ${partition}. Needs -qos wildfire"
  qosLine="#SBATCH -q wildfire"
else
  qosLine=""
fi


echo "Time Limit = ${timeLimit}"
if [ $MODE == "DNS" ] || [ $MODE == "MOVIEDNS" ] || [ $MODE == "PLOT" ] || [ $MODE == "ALL" ]; then
  echo "MODE = ${MODE}"
else
  echo "ERROR 77 - INCORRECT MODE = ${MODE}"
  exit 77
fi

OVERHEAD_THREADS=28
OVERHEAD_JOBS=14
MKL_CPUS=1

job_prefix=$(python -c "print('$job_rec'.split('/')[-1])")

sbatch_dir="log/"
sbatch_rec="${sbatch_dir}Sweep_JR${job_prefix}_${job_com// /_}"

! [[ -d "$sbatch_dir" ]] && mkdir -p "$sbatch_dir" || :

trapped() {
  echo 'TRAPPED -- QUITTING'
  exit 70
}

input_gen() {
  prefix="${1:?input_gen: PREFIX NOT PASSED}"
  restart="${2:?input_gen: RESTART NOT PASSED}"
  gamma="${3:?input_gen: GAMMA NOT PASSED}"
  eta="${4:?input_gen: ETA NOT PASSED}"
  NtsT="${5:?input_gen: NtsT NOT PASSED}"
  NT="${6:?input_gen: NT NOT PASSED}"
  Nsaves="${7:?input_gen: NSAVES NOT PASSED}"
  itseries="${8:?input_gen: ITSERIES NOT PASSED}"
  init_file="${9:?input_gen: INIT_FILE NOT PASSED}"
  iaxisym="${10:?input_gen: IAXISYM NOT PASSED}"
  ibegin="${11:?input_gen: IBEGIN NOT PASSED}"
  imode="${12:?input_gen: IMODE NOT PASSED}"
  pert="${13:?input_gen: PERT NOT PASSED}"
  out_rec="${14:?input_gen: OUT_REC NOT PASSED}"
  Bo=$(python -c "print('$prefix'.split('Bo')[-1].split('_')[0])")
  Re=$(python -c "print('$prefix'.split('Re')[-1].split('_')[0])")
  Ro=$(python -c "print('$prefix'.split('Ro')[-1].split('_')[0])")
  wf=$(python -c "print('$prefix'.split('wf')[-1].split('_')[0])")
  cat >> $out_rec << __EOF
'$prefix'     ! prefix for filenames
'$restart'    ! name of restart file
${Bo/e/d}     ! BOUSSINESQ
${Re/e/d}     ! REYNOLDS    Omega*R^2/nu
${Ro/e/d}     ! ROSSBY      Omega/Omega
${wf/e/d}     ! wf          Frequency of the sin(wot) function
${gamma/e/d}  ! GAMMA       H/R = Height/Radius (aspect ratio)
${eta/e/d}    ! eta         a/R = Knife Distance/Radius
$NtsT         ! NtsT        Number of samples for each period (used to get dt)
$NT           ! NT          Number of time periods (or time-steps, if wf=0 or Ro=0)
$Nsaves       ! Nsaves      Write Nsaves full solutions (used to calculate insec)
$itseries     ! itseries    Write in time-series-files every itseries time-steps
$init_file    ! init_file   Number of first output file
$iaxisym      ! iaxisym = m:
              !   m = 0 Axisymmetric Subspace
              !   m = 1 Full 3D Solution
              !   m > 1 m-Fourier subspace
              !   m < 0 |m|-rotoreflection
$ibegin       ! ibegin:
              !  -1 Start from solid body rotation + pert, set t=0. NOT WORKING
              !   0 Start from rest with random perturbation, set t=0.
              !   1 Continue restart solution, keep t.
              !   2 Continue restart solution, set t=0.
              !   3 Continue restart sol. with random pert., set t=0.
$imode        ! imode       Azimuthal mode to be perturbed
${pert/e/d}   ! pert        Amplitude of the vz perturbation
__EOF
  cat $out_rec
}

my_job() {
  MODE="${1:?'MODE MISSING'}"
  sourceFile="${2:?'SOURCEFILE MISSING'}"
  Bo="${3:?'BOUSSINESQ VAL MISSING'}"
  Re="${4:?'REYNOLDS VAL MISSING'}"
  Ro="${5:?'REYNOLDS VAL MISSING'}"
  wf="${6:?'FORCING FREQ VAL MISSING'}"
  gamma="${7:?'GAMMA MISSING'}"
  eta="${8:?'ETA MISSING'}"
  NtsT="${9:?'NtsT MISSING'}"
  NT="${10:?'NT MISSING'}"
  Nsaves="${11:?'NSAVES MISSING'}"
  itseries="${12:?'ITSERIES MISSING'}"
  init_file="${13:?'INIT_FILE MISSING'}"
  iaxisym="${14:?'IAXISYM MISSING'}"
  ibegin="${15:?'IBEGIN MISSING'}"
  imode="${16:?'IMODE MISSING'}"
  pert="${17:?'PERT MISSING'}"
  res_dir="${18:?'RESULT DIRECTORY MISSING'}"
  RS="${19:-NONE}"

  if [ $pert == "0e-16" ]; then
    prefix="Bo${Bo}_Re${Re}_Ro${Ro}_wf${wf}_Gamma${gamma}_eta${eta}_mode000_pert000"
  else
    prefix="Bo${Bo}_Re${Re}_Ro${Ro}_wf${wf}_Gamma${gamma}_eta${eta}_mode${imode}_pert${pert}"
  fi
  echo "$prefix"
  out_rec="${res_dir}sweep_${prefix}.out"
  ! [[ -d "$res_dir" ]] && mkdir -p "$res_dir" || :

  printf "./bin/${sourceFile}\n"
  if [ $MODE == "DNS" ] || [ $MODE == "MOVIEDNS" ] || [ $MODE == "ALL" ]; then
    printf "Computing solution for ${prefix}\n"
    ./bin/$sourceFile >> $out_rec < <(input_gen $prefix $RS $gamma\
    $eta $NtsT $NT $Nsaves $itseries $init_file $iaxisym $ibegin $imode\
    $pert $out_rec)

    if [ $MODE != "MOVIEDNS" ]; then
      mv ${prefix}*${Nsaves} "$res_dir"
      rm ${prefix}*
    fi
    mv *${prefix}* "$res_dir"
  fi
  if [ $MODE == "PLOT" ] || [ $MODE == "ALL" ]; then
    printf "Plotting timeseries for ${prefix}\n"
    ts_rec="${res_dir}ts_${prefix}"
    pycmd=$HOME/.local/opt/anaconda/bin/python
    $pycmd src/postproc/monitor.py "${ts_rec}" "${res_dir}" "${dt}"
  fi
}

export -f input_gen my_job

trap "trapped" 1 2 3 4 5 6 7 8

#SBATCH --ntasks=1
sbatch --comment="${job_com} ${job_prefix}" << EOF
#!/bin/bash
#SBATCH -p ${partition}
${qosLine}
#SBATCH --nodes=1-1
#SBATCH --exclusive
#SBATCH --mincpus=$OVERHEAD_THREADS
#SBATCH -t ${timeLimit}
#SBATCH -o "${sbatch_rec}.out"
#SBATCH -e "${sbatch_rec}.err"
#SBATCH --mail-type ALL
#SBATCH --mail-user fjcasti1@asu.edu

! [[ -d "$sbatch_dir" ]] && mkdir -p "$sbatch_dir" || :

[[ -d "../lib/" ]] && {
  export LD_LIBRARY_PATH=":$(readlink -f ../lib/):$LD_LIBRARY_PATH"
  export    LIBRARY_PATH="$LD_LIBRARY_PATH"
} || :

module load intel/2018x

export MKL_NUM_THREADS=$MKL_CPUS
ulimit -s unlimited

echo "Jobs: $job_rec"
pcmd=$HOME/.local/bin/parallel
\$pcmd -v -j $OVERHEAD_JOBS --col-sep='\s+' my_job $MODE $sourceFile :::: $job_rec
echo "FINISHED JOB"
EOF

#if [ $MODE == "PLOT" ]; then
#  if [[ $job_rec == *00 ]]; then
#    while [[ $(squeue -u fjcasti1 | wc -l) -gt "2" ]]
#    do
#      sleep 5s
#      echo "Sleeping" >> ${sbatch_rec}.out
#    done
#    python << __EOF
#import sys
#sys.path.insert(0,'/scratch/fjcasti1/generalKnifeEdge/src/PostProcessing/')
#import monitor
#monitor.collectData('dat/','*.txt','collectiveData.dat')
#__EOF
#  fi
#fi
#echo "FINISHED"


