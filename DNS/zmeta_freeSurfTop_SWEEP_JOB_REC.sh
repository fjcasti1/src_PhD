#!/usr/bin/env bash

pN=800
job_dir="job_recs/patch_${pN}/"
job_rec="${job_dir}ReSweep_JOB_REC_patch_${pN}_MASTER"

! [[ -d "$job_dir" ]] && mkdir -p "$job_dir" || :
gen_JR_line() {
# pN   Re0 ReN ReS ReE  Pe  Ca  Ro   w0  wN wS wE  gamma  NtsT NT
                                              # Nsaves dt ibegin regOpt RS
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
  gamma="${13}"
  NtsT="${14}"
  NT="${15}"
  Nsaves="${16}"
  dt="${17}"
  ibegin="${18}"
  regOpt="${19}"
  RS="${20}"
  rd=$(printf "results/runs_%s/" "$pn" )
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
gamma = "$gamma"
NtsT  = "$NtsT"
NT    = "$NT"
Nsaves= "$Nsaves"
dt    = "$dt"
ib    = "$ibegin"
regOpt= "$regOpt"
RS    = "$RS"
rd    = "$rd"

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
    print(f'{Re:s} {Pe:s} {Ca:s} {Ro:s} {wf:s} {gamma:s} {NtsT:s} {NT:s} ' +
                        f'{Nsaves:s} {dt:s} {ib:s} {regOpt:s} {rd:s} {RS:s}')
  else:
    for Frequency in range(w0,wN+1,wS):
      if Frequency < 10:
        wf = f'0{Frequency:d}{wE:s}'
      else:
        wf = f'{Frequency:d}{wE:s}'
      print(f'{Re:s} {Pe:s} {Ca:s} {Ro:s} {wf:s} {gamma:s} {NtsT:s} {NT:s} ' +
                        f'{Nsaves:s} {dt:s} {ib:s} {regOpt:s} {rd:s} {RS:s}')
__EOF
}

export -f gen_JR_line

# pN   Re0 ReN ReS ReE  Pe  Ca  Ro   w0  wN wS wE  gamma  NtsT NT
                                              # Nsaves dt ibegin regOpt RS
parallel --will-cite -j1 --col-sep='\s+' gen_JR_line :::: < <(
cat << __EOF
${pN}  1000  1000  1000 e0  0e-2 0e-2 0e-2  0 1 10 e-3  100e-2  120  100  20 5e-3  0 True NONE
__EOF
) > $job_rec
