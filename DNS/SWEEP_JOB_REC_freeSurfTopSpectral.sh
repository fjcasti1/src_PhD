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
  dt="${6:?input_gen: TIME-STEP NOT PASSED}"
  ibegin="${7:?input_gen: IBEGIN NOT PASSED}"
  regOpt="${8:?input_gen: REGOPT NOT PASSED}"
  out_rec="${9:?input_gen: OUT_REC NOT PASSED}"
  M="${10:-201}"
  N="${11:-201}"
  Re=$(python -c "print('$prefix'.split('Re')[-1].split('_')[0])")
  Pe=$(python -c "print('$prefix'.split('Pe')[-1].split('_')[0])")
  Ca=$(python -c "print('$prefix'.split('Ca')[-1].split('_')[0])")
  Ro=$(python -c "print('$prefix'.split('Ro')[-1].split('_')[0])")
  wf=$(python -c "print('$prefix'.split('wf')[-1].split('_')[0])")
  c0=$(python -c "print('$prefix'.split('co')[-1].split('_')[0])")
  gamma=$(python -c "print('$prefix'.split('Gamma')[-1].split('_')[0])")

#  npasos=$(python -c "print(int($TU/$dt))")
#  igraph=$(python -c "print(int($npasos//10))")
  (( ibegin > 0 )) && [[ -z "$restart" ]] && {
    echo 'ERROR 100 -- IBEGIN > 0 but RESTART FILE MISSING'
    exit 100
  } || :
  cat << EOF
'$prefix'        : prefix for filenames
'$restart'       : name of restart file
${Re/e/d}        : Re=Omega a^2/nu
${Pe/e/d}        : Pe
${Ca/e/d}        : Ca
${Ro/e/d}        : Ro
${wf/e/d}        : wf
${c0/e/d}        : c0, initial concentration
${gamma/e/d}     : GAMMA       H/R = Height/Radius (aspect ratio)
${M}             : Nr  = number of horizontal points
${N}             : Nz  = number of vertical points
${dt/e/d}        : dt
${NtsT/e/d}      : NtsT = Number of Time Steps per Period
${NT/e/d}        : NT   = Number of Periods
${Nsaves}        : igraph - save solution every igraph timesteps
${regOpt}        : Logical, True if regularization if top Analytical BC is deisred
1                : itseries write to ts_file every itseries timesteps
1                : init_file initial filenumber
${ibegin}        : ibegin - 0 start from rest
                          - 1 continue restart solution, keep t
                          - 2 continue restart solution, set t=0
EOF
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
  NtsT="${10:?'NtsT MISSING'}"
  NT="${11:?'NT MISSING'}"
  Nsaves="${12:?'N SAVES MISSING'}"
  dt="${13:?'TIME-STEP VAL MISSING'}"
  ibegin="${14:?'IBEGIN (STATE) MISSING'}"
  regOpt="${15:?'REG OPT MISSING'}"
  res_dir="${16:?'RESULT DIRECTORY MISSING'}"
  restart="${17:-'NONE'}"
  M="${18:-201}"
  N="${19:-201}"

  prefix="Re${Re}_Pe${Pe}_Ca${Ca}_Ro${Ro}_wf${wf}_co${c0}_Gamma${gamma}_NtsT${NtsT}_NT${NT}"
  out_rec="${res_dir}sweep_${prefix}.out"
  ! [[ -d "$res_dir" ]] && mkdir -p "$res_dir" || :

  printf "Running binary file ${sourceFile}\n"
  if [ $MODE == "DNS" ] || [ $MODE == "MOVIEDNS" ] || [ $MODE == "ALL" ]; then
    printf "Computing solution for ${prefix}\n"
    $sourceFile > $out_rec < <(
      input_gen $prefix $restart $NtsT $NT $Nsaves $dt $ibegin $regOpt $out_rec $M $N)

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
#import monitor_freeSurfTop
#monitor_freeSurfTop.collectData('dat/','*.txt','collectiveData.dat')
#__EOF
#  fi
#fi
#echo "FINISHED"


