#!/usr/bin/env bash
job_rec="${1:?'JOB RECORD NOT SUPPLIED'}"
sourceFile="${2:?'SOURCE FILE NOT SUPPLIED'}"
timeLimit="${3:?'TIME LIMIT NOT SUPPLIED'}"
partition="${4:?'PARTITION NOT SUPPLIED'}"
MODE="${5:?'MODE NOT SUPPLIED'}"
job_com="${6:-main}"

if [ -x ${sourceFile} ]; then
  echo "Running binary file ${sourceFile}"
else
  echo "ERROR 66 - INCORRECT BINARY FILE ${sourceFile}"
  exit 66
fi

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
  NtsT="${3:?input_gen: NtsT NOT PASSED}"
  NT="${4:?input_gen: NT NOT PASSED}"
  Nsaves="${5:?input_gen: NSAVES NOT PASSED}"
  itseries="${6:?input_gen: ITSERIES NOT PASSED}"
  init_file="${7:?input_gen: INIT_FILE NOT PASSED}"
  iaxisym="${8:?input_gen: IAXISYM NOT PASSED}"
  ibegin="${9:?input_gen: IBEGIN NOT PASSED}"
  imode="${10:?input_gen: IMODE NOT PASSED}"
  pert="${11:?input_gen: PERT NOT PASSED}"
  out_rec="${12:?input_gen: OUT_REC NOT PASSED}"
  dt="${13:?input_gen: TIME-STEP NOT PASSED}"
  Re=$(python -c "print('$prefix'.split('Re')[-1].split('_')[0])")
  Pe=$(python -c "print('$prefix'.split('Pe')[-1].split('_')[0])")
  Ca=$(python -c "print('$prefix'.split('Ca')[-1].split('_')[0])")
  Ro=$(python -c "print('$prefix'.split('Ro')[-1].split('_')[0])")
  wf=$(python -c "print('$prefix'.split('wf')[-1].split('_')[0])")
  c0=$(python -c "print('$prefix'.split('co')[-1].split('_')[0])")
  gamma=$(python -c "print('$prefix'.split('Gamma')[-1].split('_')[0])")
  (( ibegin > 0 )) && [[ -z "$restart" ]] && {
    echo 'ERROR 100 -- IBEGIN > 0 but RESTART FILE MISSING'
    exit 100
  } || :

  cat >> $out_rec << __EOF
'$prefix'     ! prefix for filenames
'$restart'    ! name of restart file
${Re/e/d}     ! REYNOLDS    Omega*R^2/nu
${Pe/e/d}     ! PECLET      Omega*R^2/D_s
${Ca/e/d}     ! CAPILLARY   mu*Omega*R/sigma_0
${Ro/e/d}     ! ROSSBY      Amplitud/Omega
${wf/e/d}     ! wf          Frequency of the sin(wf*t) function
${c0/e/d}     ! c0          Initial concentration
${gamma/e/d}  ! GAMMA       H/R = Height/Radius (aspect ratio)
$dt           ! dt          time-step
$NtsT         ! NtsT        Number of samples for each forcing period (used to get dt)
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
  Re="${3:?'REYNOLDS VAL MISSING'}"
  Pe="${4:?'PECLET VAL MISSING'}"
  Ca="${5:?'CAPILARY VAL MISSING'}"
  Ro="${6:?'ROSSBY VAL MISSING'}"
  wf="${7:?'FORCING FREQ VAL MISSING'}"
  c0="${8:?'INITIAL CONCENTRATION MISSING'}"
  gamma="${9:?'GAMMA MISSING'}"
  dt="${10:?'TIME-STEP VAL MISSING'}"
  NtsT="${11:?'NtsT MISSING'}"
  NT="${12:?'NT MISSING'}"
  Nsaves="${13:?'NSAVES MISSING'}"
  itseries="${14:?'ITSERIES MISSING'}"
  init_file="${15:?'INIT_FILE MISSING'}"
  iaxisym="${16:?'IAXISYM MISSING'}"
  ibegin="${17:?'IBEGIN MISSING'}"
  imode="${18:?'IMODE MISSING'}"
  pert="${19:?'PERT MISSING'}"
  res_dir="${20:?'RESULT DIRECTORY MISSING'}"
  restart="${21:-'NONE'}"

  if [ $pert == "0e-16" ]; then
    prefix="Re${Re}_Pe${Pe}_Ca${Ca}_Ro${Ro}_wf${wf}_co${c0}_Gamma${gamma}_mode000_pert000"
  else
    prefix="Re${Re}_Pe${Pe}_Ca${Ca}_Ro${Ro}_wf${wf}_co${c0}_Gamma${gamma}_mode${imode}_pert${pert}"
  fi
  echo "$prefix"
  out_rec="${res_dir}sweep_${prefix}.out"
  ! [[ -d "$res_dir" ]] && mkdir -p "$res_dir" || :

  printf "Running binary file ${sourceFile}\n"
  if [ $MODE == "DNS" ] || [ $MODE == "MOVIEDNS" ] || [ $MODE == "ALL" ]; then
    printf "Computing solution for ${prefix}\n"
    $sourceFile >> $out_rec < <(input_gen $prefix $restart $NtsT $NT $Nsaves\
      $itseries $init_file $iaxisym $ibegin $imode $pert $out_rec $dt)

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
    $pycmd src_PhD/postproc/monitor_freeSurfTop.py "${ts_rec}" "${dt}"
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


