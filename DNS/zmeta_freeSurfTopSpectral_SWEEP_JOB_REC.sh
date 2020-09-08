#!/usr/bin/env bash
pN=000
job_dir="job_recs/patch_${pN}/"
job_rec="${job_dir}ReSweep_JOB_REC_patch_${pN}_MASTER"

! [[ -d "$job_dir" ]] && mkdir -p "$job_dir" || :
gen_JR_line() {
  pn="$1"
  Re0="$2"
  ReN="$3"
  ReS="$4"
  ReE="$5"
  Pe="$6"
  Ca="$7"
  Ro="$8"
  w0="$9"
  wN="${10}"
  wS="${11}"
  wE="${12}"
  c0="${13}"
  gamma="${14}"
  dt="${15}"
  NtsT="${16}"
  NT="${17}"
  Nsaves="${18}"
  itseries="${18}"
  init_file="${19}"
  iaxisym="${20}"
  ibegin="${21}"
  imode="${22}"
  pert="${23}"
  RS="${24}"
  res_dir=$(printf "results/runs_%s/" "$pn")
  python << __EOF
from numpy import pi
Re0   =  $Re0
ReN   =  $ReN
ReS   =  $ReS
ReE   = "$ReE"
Pe    = "$Pe"
Ca    = "$Ca"
Ro    = "$Ro"
w0    =  $w0
wN    =  $wN
wS    =  $wS
wE    = "$wE"
c0    = "$c0"
gamma = "$gamma"
dt    = "$dt"
NtsT  = "$NtsT"
NT    = "$NT"
Nsaves= "$Nsaves"
itseries = "$itseries"
init_file = "$init_file"
iaxisym = "$iaxisym"
ibegin = "$ibegin"
imode = "$imode"
pert = "$pert"
RS = "$RS"
res_dir = "$res_dir"

RS    = "" if RS == 'NONE' else RS

if ib == 0 and RS != '':
  raise ValueError("ibegin is 0 but restart is specified!")
elif ib == 1 and RS != 'RS':
  raise ValueError("ibegin is 1 but restart is not specified!")
elif ib == 2 and RS != 'RS':
  raise ValueError("ibegin is 2 but restart is not specified!")

for Reynolds in range(Re0,ReN+1,ReS):
  if ReE == 'e0':
    Re = f'{Reynolds:d}'
  else:
    Re = f'{Reynolds:d}{ReE:s}'
  if w0 == wN:
    Frequency = w0
    Frequency = Frequency*float('1'+wE)
    wf = str(int(Frequency*1e18))+'e-18'
    print(f'{Re:s} {Pe:s} {Ca:s} {Ro:s} {wf:s} {c0:s} {gamma:s} ' +
          f'{dt:s} {NtsT:s} {NT:s} {Nsaves:s} {itseries:s} {init_file:s} ' +
          f'{iaxisym:s} {ibegin:s} {imode:s} {pert:s} {res_dir:s} {RS:s}')
  else:
    for Frequency in range(w0,wN+1,wS):
      if Frequency < 10:
        wf = f'0{Frequency:d}{wE:s}'
      else:
        wf = f'{Frequency:d}{wE:s}'
      print(f'{Re:s} {Pe:s} {Ca:s} {Ro:s} {wf:s} {c0:s} {gamma:s} ' +
            f'{dt:s} {NtsT:s} {NT:s} {Nsaves:s} {itseries:s} {init_file:s} ' +
            f'{iaxisym:s} {ibegin:s} {imode:s} {pert:s} {res_dir:s} {RS:s}')
__EOF
}

export -f gen_JR_line

# pN   Re0 ReN ReS ReE  Pe  Ca  Ro   w0  wN wS wE  gamma  NtsT NT
                                              # Nsaves dt ibegin regOpt RS
parallel --will-cite -j1 --col-sep='\s+' gen_JR_line :::: < <(
cat << __EOF
__EOF
) > $job_rec