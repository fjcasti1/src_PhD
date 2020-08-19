#!/usr/bin/env python
import sys, os
import numpy as np
from scipy.signal import find_peaks
from math import ceil
from pylab import detrend,fft,savefig
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from scipy.signal import blackman as blk
from glob import glob
import pandas as pd


def fig_dir(f):
  #a = f.split('/')[0]
  #return f'fig/{a:s}/'
  return f'fig/'

def long_basename(f):
  return f.split('/')[-1]

def ts_basename(f):
  return f.split('_NtsT')[0]

def fft_basename(f):
  return f.replace('ts_','fft_')

def orbit_basename(f):
  return f.replace('ts_','orbit_')

def parse_token(bn,token):
  return bn.split(token)[-1].split('_')[0]

def pktopkAmp(S,M=0,F=0.9):
  """
    Help here
  """
  if M == 0:
    M = len(S)
  thresUp  = np.mean(S)+F*(np.max(S)-np.mean(S))
  thresDwn = np.mean(S)-F*(np.mean(S)-np.min(S))
  peaks, _   = find_peaks(S[-M:], height=thresUp)
  valleys, _ = find_peaks(-S[-M:],height=-thresDwn)
  pkpkAmp = np.mean(S[-M:][peaks]) - np.mean(S[-M:][valleys])
  relErr = ((np.std(S[-M:][peaks])**2 + np.std(S[-M:][valleys])**2)**0.5)/pkpkAmp
  return (pkpkAmp, relErr)

def findIndex(df,Bo,Re,alpha,wf):
  cond = ( (df['Re']==Re) & (df['Bo']==Bo) & (df['alpha']==alpha) &
      (df['w_f']==wf) )
  return df.index[cond].tolist()

def addRow(df,runs,Bo,Re,alpha,wf,NtsT,NT,wFFT,AEk,AvgEk):
  df = df.append({'runs_#':runs, 'Bo':Bo, 'Re':Re, 'alpha':alpha, 'w_f':wf,
    'NtsT':NtsT, 'NT':NT, 'w*':wFFT, 'stdEk': AEk, 'AvgEk': AvgEk}, ignore_index=True)
  return df

def replaceRow(df,index,runs,Bo,Re,alpha,wf,NtsT,NT,wFFT,stdEk,AvgEk):
  df.loc[index,['runs_#','Bo','Re','alpha','w_f','NtsT','NT','w*',
    'stdEk','AvgEk']]=[runs,Bo,Re,alpha,wf,NtsT,NT,wFFT,stdEk,AvgEk]
  return None

def collectData(DAT_DIR,infiles,outfile):
  df = pd.DataFrame(columns=['runs_#','Bo','Re','alpha','w_f','NtsT',
    'NT','w*','stdEk','AvgEk'])
  if os.path.exists(DAT_DIR+outfile):
    df = pd.read_csv(DAT_DIR+outfile, sep=' ', dtype=object)
  for infile in glob(DAT_DIR+infiles):
    with open(infile,'r') as f:
      f.readline()
      try:
        params = f.readline().strip('\n').split()
        runs      = params[0]
        Bo        = params[1]
        Re        = params[2]
        alpha     = params[3]
        wf        = params[4]
        NtsT      = params[5]
        NT        = params[6]
        wFFT      = params[7]
        stdEk     = params[8]
        AvgEk     = params[9]
      except Exception as ex:
        print('Exception reading line: ', ex)
      filterIndex = findIndex(df,Bo,Re,alpha,wf)
      if filterIndex and runs >= df.loc[filterIndex,'runs_#'].values:
        replaceRow(df,filterIndex,runs,Bo,Re,alpha,wf,NtsT,NT,wFFT,stdEk,AvgEk)
      elif not filterIndex:
        df = addRow(df,runs,Bo,Re,alpha,wf,NtsT,NT,wFFT,stdEk,AvgEk)
      f.close()
    os.remove(infile)

  with open(DAT_DIR+outfile,'w') as outfile:
    df.to_csv(outfile,header=True,index=False,sep=' ')
    outfile.close()
  return None

def main():
  f       = sys.argv[1]
  res_dir = sys.argv[2]
  dtinput = sys.argv[3]

  FIG_DIR = fig_dir(res_dir)
  DAT_DIR = f'dat/'

  longbn  = long_basename(f)
  tsbn    = ts_basename(longbn)
  fftbn   = fft_basename(tsbn)
  orbitbn = orbit_basename(tsbn)
  tokens  = ['Re','Bo','alpha','wf','NtsT','NT']
  try:
    values = [ parse_token(longbn,token) for token in tokens ]
    Re   = values[0]
    Bo   = values[1]
    alpha= values[2]
    wf   = values[3]  # Forcing Angular Freq
    NtsT = int(values[4])
    NT   = int(values[5])
    runs = f.split('runs_')[-1].split('/')[0]
  except Exception as ex:
    print('Exception in parse token: ', ex)
  if float(wf)>0:
    Period = 2*np.pi/float(wf)
    dt     = Period/NtsT
    Nsteps = float(NT*NtsT)
  else:
    print('wf not greater than 0. wf = ', wf) # COMPLETE THIS!

#  if TU<0:
#    Nsteps = -TU
#    if wf!=0:
#      dt     = float(1/(wf*Nsteps))
#    elif wf==0:
#      dt = float(dtinput)
#  else:
#    dt = float(dtinput)
#    Nsteps = float(TU/dt)

  title_string = longbn.replace('_',' ')
  t,Ek,Eg,Ew,ur,uw,uz = np.loadtxt(f).T
  os.makedirs(FIG_DIR,exist_ok=True)


#############
# Plot ffts #
#############
  P = 100
  M = int(len(Ek)*P/100)
  T = M*dt      # Period ??
  w0  = 2*np.pi/T  # Natural Frequency??

  AEk = Ek[-M:].std() # Amplitud of Oscillation
  fftEk  = abs(fft(detrend(Ek[-M:])*blk(M))[:M//2]) # FFT with Blackman filter [array]
  wMEk = w0*fftEk.argmax() # Compute dominant frequency

  AEg = Eg[-M:].std() # Amplitud of Oscillation
  fftEg  = abs(fft(detrend(Eg[-M:])*blk(M))[:M//2]) # FFT with Blackman filter [array]
  wMEg = w0*fftEg.argmax() # Compute dominant frequency

  AEw = Ew[-M:].std() # Amplitud of Oscillation
  fftEw  = abs(fft(detrend(Ew[-M:])*blk(M))[:M//2]) # FFT with Blackman filter [array]
  wMEw = w0*fftEw.argmax() # Compute dominant frequency

  Aur = ur[-M:].std() # Amplitud of Oscillation
  fftur  = abs(fft(detrend(ur[-M:])*blk(M))[:M//2]) # FFT with Blackman filter [array]
  wMur = w0*fftur.argmax() # Compute dominant frequency

  Auw = uw[-M:].std() # Amplitud of Oscillation
  fftuw  = abs(fft(detrend(uw[-M:])*blk(M))[:M//2]) # FFT with Blackman filter [array]
  wMuw = w0*fftuw.argmax() # Compute dominant frequency

  Auz = uz[-M:].std() # Amplitud of Oscillation
  fftuz  = abs(fft(detrend(uz[-M:])*blk(M))[:M//2]) # FFT with Blackman filter [array]
  wMuz = w0*fftuz.argmax() # Compute dominant frequency

  wFFT = min([wMEk,wMEg,wMEw,wMur,wMuw,wMuz])

  wLim = 2
  AnotationSize = 15
  xPosText = 0.25
  yPosText = 0.92
  ticksize = 12
  labelsize = 18
  labelpadx = 3
  labelpady = 16

  fig, axes = plt.subplots(nrows=2,ncols=3,figsize=(14,9)) # Create canvas & axes
  ## Global Kinetic Energy FFT
  axes[0,0].semilogy(w0*np.arange(len(fftEk)),fftEk,'k-')
  axes[0,0].annotate('$\omega^*$ = {:f}'.format(wMEk), xy=(wMEk, fftEk.max()),
          xycoords='data', xytext=(xPosText,yPosText), textcoords='axes fraction',
          size=AnotationSize, arrowprops=dict(arrowstyle="->"))
  axes[0,0].set_xlabel('$\omega$',fontsize=labelsize,labelpad=labelpadx)
  axes[0,0].set_ylabel('$|\hat{E}_k|$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  axes[0,0].set_xlim(0,wLim)
  axes[0,0].tick_params(labelsize=ticksize)
  ## Global Angular Momentum FFT
  axes[0,1].semilogy(w0*np.arange(len(fftEw)),fftEw,'k-')
  axes[0,1].annotate('$\omega^*$ = {:f}'.format(wMEw), xy=(wMEw, fftEw.max()),
          xycoords='data', xytext=(xPosText,yPosText), textcoords='axes fraction',
          size=AnotationSize, arrowprops=dict(arrowstyle="->"))
  axes[0,1].set_xlabel('$\omega$',fontsize=labelsize,labelpad=labelpadx)
  axes[0,1].set_ylabel('$|\hat{E}_w|$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  axes[0,1].set_xlim(0,wLim)
  axes[0,1].tick_params(labelsize=ticksize)
  ## Global Enstrophy FFT
  axes[0,2].semilogy(w0*np.arange(len(fftEg)),fftEg,'k-')
  axes[0,2].annotate('$\omega^*$ = {:f}'.format(wMEg), xy=(wMEg, fftEk.max()),
          xycoords='data', xytext=(xPosText,yPosText), textcoords='axes fraction',
          size=AnotationSize, arrowprops=dict(arrowstyle="->"))
  axes[0,2].set_xlabel('$\omega$',fontsize=labelsize,labelpad=labelpadx)
  axes[0,2].set_ylabel('$|\hat{E}_{\gamma}|$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  axes[0,2].set_xlim(0,wLim)
  axes[0,2].tick_params(labelsize=ticksize)
  ## Local Radial Velocity FFT
  axes[1,0].semilogy(w0*np.arange(len(fftur)),fftur,'k-')
  axes[1,0].annotate('$\omega^*$ = {:f}'.format(wMur), xy=(wMur, fftur.max()),
          xycoords='data', xytext=(xPosText,yPosText), textcoords='axes fraction',
          size=AnotationSize, arrowprops=dict(arrowstyle="->"))
  axes[1,0].set_xlabel('$\omega$',fontsize=labelsize,labelpad=labelpadx)
  axes[1,0].set_ylabel('$|\hat{u}_r|$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  axes[1,0].set_xlim(0,wLim)
  axes[1,0].tick_params(labelsize=ticksize)
  ## Local Azimuthal Velocity FFT
  axes[1,1].semilogy(w0*np.arange(len(fftuw)),fftuw,'k-')
  axes[1,1].annotate('$\omega^*$ = {:f}'.format(wMuw), xy=(wMuw, fftuw.max()),
          xycoords='data', xytext=(xPosText,yPosText), textcoords='axes fraction',
          size=AnotationSize, arrowprops=dict(arrowstyle="->"))
  axes[1,1].set_xlabel('$\omega$',fontsize=labelsize,labelpad=labelpadx)
  axes[1,1].set_ylabel(r'$|\hat{u}_{\theta}|$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  axes[1,1].set_xlim(0,wLim)
  axes[1,1].tick_params(labelsize=ticksize)
  ## Local Axial Velocity FFT
  axes[1,2].semilogy(w0*np.arange(len(fftuz)),fftuz,'k-')
  axes[1,2].annotate('$\omega^*$ = {:f}'.format(wMuz), xy=(wMuz, fftuz.max()),
          xycoords='data', xytext=(xPosText,yPosText), textcoords='axes fraction',
          size=AnotationSize, arrowprops=dict(arrowstyle="->"))
  axes[1,2].set_xlabel('$\omega$',fontsize=labelsize,labelpad=labelpadx)
  axes[1,2].set_ylabel('$|\hat{u}_z|$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  axes[1,2].set_xlim(0,wLim)
  axes[1,2].tick_params(labelsize=ticksize)

  fig.tight_layout()
  savefig(f'{FIG_DIR:s}{fftbn:s}.png')
  plt.close()


####################
# Plot time series #
####################
  P = 5 # last P% of the time series
  Nperiods = 4
  if float(wf)>0:
    M = int(Nperiods*NtsT)
#  else:                               COMPLETE THIS PART
#    if wFFT == 0:
#      M = int(len(t)*P/100)
#    else:
#      TUmin = Nperiods*2*np.pi/wFFT
#      M = ceil(TUmin/dt)
  ticksize = 12
  labelsize = 18
  labelpadx = 3
  labelpady = 10
  w = 1+float(alpha)*np.cos(float(wf)*t[-M:])

  fig, axes = plt.subplots(nrows=2,ncols=3,figsize=(14,9)) # Create canvas & axes
  ## Global Kinetic Energy Time Series
  axes[0,0].plot(t[-M:]/Period,Ek[-M:],'r-')
  axes[0,0].set_xlabel('$t/T_f$',fontsize=labelsize,labelpad=labelpadx)
  axes[0,0].set_ylabel('$E_k$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  axes[0,0].tick_params(labelsize=ticksize)
  ax2 = axes[0,0].twinx()
  ax2.plot(t[-M:]/Period,w, color='tab:blue')
  ax2.set_ylabel('$\omega$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  ax2.tick_params(labelsize=ticksize)
  ## Global Angular Momentum Time Series
  axes[0,1].plot(t[-M:]/Period,Ew[-M:],'r-')
  axes[0,1].set_xlabel('$t/T_f$',fontsize=labelsize,labelpad=labelpadx)
  axes[0,1].set_ylabel('$E_{\omega}$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  axes[0,1].tick_params(labelsize=ticksize)
  ax2 = axes[0,1].twinx()
  ax2.plot(t[-M:]/Period,w, color='tab:blue')
  ax2.set_ylabel('$\omega$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  ax2.tick_params(labelsize=ticksize)
  ## Global Enstrophy Time Series
  axes[0,2].plot(t[-M:]/Period,Eg[-M:],'r-')
  axes[0,2].set_xlabel('$t/T_f$',fontsize=labelsize,labelpad=labelpadx)
  axes[0,2].set_ylabel('$E_{\gamma}$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  axes[0,2].tick_params(labelsize=ticksize)
  ax2 = axes[0,2].twinx()
  ax2.plot(t[-M:]/Period,w, color='tab:blue')
  ax2.set_ylabel('$\omega$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  ax2.tick_params(labelsize=ticksize)
  ## Local Radial Velocity Time Series
  axes[1,0].plot(t[-M:]/Period,ur[-M:],'r-')
  axes[1,0].set_xlabel('$t/T_f$',fontsize=labelsize,labelpad=labelpadx)
  axes[1,0].set_ylabel('$u_r$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  axes[1,0].tick_params(labelsize=ticksize)
  ax2 = axes[1,0].twinx()
  ax2.plot(t[-M:]/Period,w, color='tab:blue')
  ax2.set_ylabel('$\omega$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  ax2.tick_params(labelsize=ticksize)
  ## Local Azimuthal Velocity Time Series
  axes[1,1].plot(t[-M:]/Period,uw[-M:],'r-')
  axes[1,1].set_xlabel('$t/T_f$',fontsize=labelsize,labelpad=labelpadx)
  axes[1,1].set_ylabel(r'$u_{\theta}$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  axes[1,1].tick_params(labelsize=ticksize)
  ax2 = axes[1,1].twinx()
  ax2.plot(t[-M:]/Period,w, color='tab:blue')
  ax2.set_ylabel('$\omega$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  ax2.tick_params(labelsize=ticksize)
  ## Local Axial Velocity Time Series
  axes[1,2].plot(t[-M:]/Period,uz[-M:],'r-')
  axes[1,2].set_xlabel('$t/T_f$',fontsize=labelsize,labelpad=labelpadx)
  axes[1,2].set_ylabel('$u_z$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  axes[1,2].tick_params(labelsize=ticksize)
  ax2 = axes[1,2].twinx()
  ax2.plot(t[-M:]/Period,w, color='tab:blue')
  ax2.set_ylabel('$\omega$',rotation=0,fontsize=labelsize,labelpad=labelpady)
  ax2.tick_params(labelsize=ticksize)

  fig.tight_layout()
  fig.savefig(f'{FIG_DIR:s}{tsbn:s}.png')
  plt.close()

#####################
# Plot phase orbits #
#####################
  elevation = 30
  theta0 = 10
  dtheta = 35
  ticksize = 0.1
  labelsize = 16
  labelpadx = 0
  labelpady = 0
  labelpadz = 0

  fig = plt.figure(figsize=(14,10)) # Create canvas
  ## Plot Global Orbit
  for j in range(1,4):
      ax = fig.add_subplot(2,3,j,projection='3d')
      ax.xaxis.set_rotate_label(False)  # disable automatic rotation
      ax.yaxis.set_rotate_label(False)  # disable automatic rotation
      ax.zaxis.set_rotate_label(False)  # disable automatic rotation
      ax.plot(Eg,Ew,Ek,'g-')
      ax.set_xlabel('$E_{\gamma}$',fontsize=labelsize,labelpad=labelpadx)
      ax.set_ylabel('$E_{\omega}$',rotation=0,fontsize=labelsize,labelpad=labelpady)
      ax.set_zlabel('$E_k$',rotation=0,fontsize=labelsize,labelpad=labelpadz)
      ax.tick_params(labelsize=ticksize)
      ax.view_init(elevation,theta0+(j-1)*dtheta)
  ## Plot Local Orbit
  for j in range(1,4):
      ax = fig.add_subplot(2,3,j+3,projection='3d')
      ax.xaxis.set_rotate_label(False)  # disable automatic rotation
      ax.yaxis.set_rotate_label(False)  # disable automatic rotation
      ax.zaxis.set_rotate_label(False)  # disable automatic rotation
      ax.plot(ur,uw,uz,'b-')
      ax.set_xlabel('$u_r$',fontsize=labelsize,labelpad=labelpadx)
      ax.set_ylabel(r'$u_{\theta}$',rotation=0,fontsize=labelsize,labelpad=labelpady)
      ax.set_zlabel('$u_z$',rotation=0,fontsize=labelsize,labelpad=labelpadz)
      ax.tick_params(labelsize=ticksize)
      ax.view_init(elevation,theta0+(j-1)*dtheta)

  fig.tight_layout()
  savefig(f'{FIG_DIR:s}{orbitbn:s}.png')
  plt.close()


##############
# Write Data #
##############
  #(pkpkAmpEk, relErrEk) = pktopkAmp(Ek)
  AvgEk = sum(Ek[-NtsT:])*dt/Period
  dataFile = longbn+'.txt'
  df = pd.DataFrame(columns=['runs_#','Bo','Re','alpha','w_f','NtsT','NT','w*','stdEk','AvgEk'])
  df = addRow(df,runs,Bo,Re,alpha,wf,NtsT,NT,wFFT,AEk,AvgEk)

  with open(os.path.join(DAT_DIR, dataFile),'w') as outfile:
    df.to_csv(outfile,header=True,index=False,sep=' ')
    outfile.close()

  return None


if __name__ == '__main__':
  main()

