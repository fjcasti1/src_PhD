#!/usr/bin/env bash
pN=900
job_dir="job_recs/patch_${pN}/"
job_rec="${job_dir}ReSweep_JOB_REC_patch_${pN}_MASTER"

! [[ -d "$job_dir" ]] && mkdir -p "$job_dir" || :
gen_JR_line() {
  pn="$1"
  Bo="$2"
  Re0="$3"
  ReN="$4"
  ReS="$5"
  ReE="$6"
  Ro="$7"
  w0="$8"
  wN="$9"
  wS="${10}"
  wE="${11}"
  gamma="${12}"
  eta="${13}"
  NtsT="${14}"
  NT="${15}"
  Nsaves="${16}"
  itseries="${17}"
  init_file="${18}"
  iaxisym="${19}"
  ibegin="${20}"
  imode="${21}"
  pert="${22}"
  RS="${23}"
  res_dir=$(printf "results/runs_%s/" "$pn")
  python << __EOF
from numpy import pi
Bo        = "$Bo"
Re0       =  $Re0
ReN       =  $ReN
ReS       =  $ReS
ReE       = "$ReE"
Ro        = "$Ro"
w0        =  $w0
wN        =  $wN
wS        =  $wS
wE        = "$wE"
gamma     = "$gamma"
eta       = "$eta"
NtsT      = "$NtsT"
NT        = "$NT"
Nsaves    = "$Nsaves"
itseries  = "$itseries"
init_file = "$init_file"
iaxisym   = "$iaxisym"
ib        = "$ibegin"
imode     = "$imode"
pert      = "$pert"
RS        = "$RS"
res_dir   = "$res_dir"

RS    = "" if RS == 'NONE' else RS

if ib == -1 and RS != '':
  raise ValueError("ibegin is -1 but restart is specified!")
elif ib == 0 and RS != '':
  raise ValueError("ibegin is 0 but restart is specified!")
elif ib == 1 and RS != 'RS':
  raise ValueError("ibegin is 1 but restart is not specified!")
elif ib == 2 and RS != 'RS':
  raise ValueError("ibegin is 2 but restart is not specified!")
elif ib == 3 and RS != 'RS':
  raise ValueError("ibegin is 3 but restart is not specified!")

for Reynolds in range(Re0,ReN+1,ReS):
  if ReE == 'e0':
    Re = f'{Reynolds:d}'
  else:
    Re = f'{Reynolds:d}{ReE:s}'
  if w0 == wN:
    Frequency = w0
    Frequency = Frequency*float('1'+wE)
    wf = str(int(Frequency*1e18))+'e-18'
    print(f'{Bo:s} {Re:s} {Ro:s} {wf:s} {gamma:s} {eta:s} {NtsT:s} {NT:s} ' +
      f'{Nsaves:s} {itseries:s} {init_file:s} {iaxisym:s} {ib:s} {imode:s} ' +
      f'{pert:s} {res_dir:s} {RS:s}')
  else:
    for Frequency in range(w0,wN+1,wS):
      if Frequency < 10:
        wf = f'0{Frequency:d}{wE:s}'
      else:
        wf = f'{Frequency:d}{wE:s}'
      print(f'{Bo:s} {Re:s} {Ro:s} {wf:s} {gamma:s} {eta:s} {NtsT:s} {NT:s} ' +
        f'{Nsaves:s} {itseries:s} {init_file:s} {iaxisym:s} {ib:s} {imode:s} ' +
        f'{pert:s} {res_dir:s} {RS:s}')
__EOF
}

export -f gen_JR_line

# pN   Bo  Re0 ReN ReS ReE  Ro   w0  wN wS wE  gamma eta  NtsT NT Nsaves
                            # itseries init_file iaxisym ibegin imode pert RS
parallel --will-cite -j1 --col-sep='\s+' gen_JR_line :::: < <(
cat << __EOF
${pN}  0e0   7000  9000  1000 e0  0e-2  0 1 10 e-3  100e-2  50e-2  10000 200 20 10 1 0 0 0 0e-16
${pN}  5e2   7000  9000  1000 e0  0e-2  0 1 10 e-3  100e-2  50e-2  10000 200 20 10 1 0 0 0 0e-16
${pN}  5e1   7000  9000  1000 e0  0e-2  0 1 10 e-3  100e-2  50e-2  10000 200 20 10 1 0 0 0 0e-16
__EOF
) > $job_rec
#${pN}  5e1  10000 14000 1000 e0  0e-2  0 1 10 e-3  075e-2  50e-2  10000 200 20 10 1 0 2 0 0e-16
#${pN}  5e1   6000  6500  100 e0  0e-2  0 1 10 e-3  085e-2  50e-2  10000 200 20 10 1 0 2 0 0e-16
#${pN}  5e1   5200  5700  100 e0  0e-2  0 1 10 e-3  095e-2  50e-2  10000 200 20 10 1 0 2 0 0e-16
#${pN}  5e1   5000  5700  100 e0  0e-2  0 1 10 e-3  105e-2  50e-2  10000 200 20 10 1 0 2 0 0e-16
#${pN}  5e1   6000  6500  100 e0  0e-2  0 1 10 e-3  115e-2  50e-2  10000 200 20 10 1 0 2 0 0e-16
#${pN}  5e1   7000  7500  100 e0  0e-2  0 1 10 e-3  125e-2  50e-2  10000 200 20 10 1 0 2 0 0e-16
