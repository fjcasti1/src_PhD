import glob,natsort
import multiprocessing as mp
import pandas as pd
import sys,os
from numpy import *
from pylab import *
from scipy.signal import blackman as blk

NPROCS = 10
PLOT_THRESHOLD = 1e-25
LEGEND_THRESHOLD = 1e-30
titlesize = 18
legendfontsize = 12
labelsize = 20
labelpadx = 6
labelpady = 25
colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple',
          'tab:brown','tab:pink','tab:cyan','gold','lime','teal','violet',
          'chocolate','tomato','olive','magenta'] # length = 16

def get_figdir(f):
    fig_dir = f'fig/'
    return fig_dir

def get_basename(f):
    return f.split('/ts')[-1]

def get_token(bn,token):
    return bn.split(token)[1].split('_')[0]

def parse_header(header):
    Nrows   = int(header.split('Nrows=')[-1].split(',')[0].strip())
    Ncols   = int(header.split('Ncols=')[-1].split(',')[0].strip())
    Nmodes  = int(header.split('Nmodes=')[-1].split(',')[0].strip())
    Nprobes = int(header.split('Nprobes=')[-1].split(',')[0].strip())
    return Nrows,Ncols,Nmodes,Nprobes

def get_figname(f):
    bn = get_basename(f)
    figname = 'fig/{:s}.png'.format(bn)
    return figname

def unpack_time_series(d,Ncols,Nmodes,Nprobes):
    D = array(d).reshape(len(d)//Ncols,Ncols).T
    t = D[0,:]
    Ekm   = zeros((len(t),Nmodes ))
    phspt = zeros((len(t),Nprobes))
    for i in range(0,Nmodes): #range already goes to Nmodes-1
        Ekm[:,i] = D[i+1,:]
    for i in range(0,Nprobes):#range already goes to Nprobes-1
        phspt[:,i] = D[Nmodes+i+1,:]
    return t,Ekm,phspt

def get_fields(f):
    with open(f,'rb') as ff:
        header = str(ff.read(256))
    Nrows,Ncols,Nmodes,Nprobes = parse_header(header)
    d = memmap(f,dtype=double,offset=256,mode='r')
    t,Ekm,phspt = unpack_time_series(d,Ncols,Nmodes,Nprobes)
    del d
    basename = get_basename(f)
    Bo    = get_token(basename,'Bo')
    Re    = get_token(basename,'Re')
    Ro    = get_token(basename,'Ro')
    wf    = get_token(basename,'wf')
    gamma = get_token(basename,'Gamma')
#  # Save fields in dictionary
    fields = {
        'f'     : f,
        'bn'    : basename,
        'Bo'    : Bo,
        'Re'    : Re,
        'Ro'    : Ro,
        'wf'    : wf,
        'gamma' : gamma,
        't'     : t,
        'Ekm'   : Ekm,
        'phspt' : phspt
    }
    print('parsed: {:s}'.format(f))
    return fields

def plotter(fields,fig_dir):
    bn    = fields['bn']
    Bo    = float(fields['Bo'])
    Re    = float(fields['Re'])
    Ro    = float(fields['Ro'])
    wf    = float(fields['wf'])
    gamma = float(fields['gamma'])
    t     = fields['t']
    Ekm   = fields['Ekm']
    # Calculate Total Kinetic Energy
    Ek    = (sum(Ekm,axis=1)) # Should be multiplied by 2pi?
    # Calculate Boundary Condition of the Knife
    w  = 1+Ro*sin(wf*t)

##################
# Kinetic Energy #
##################
    P = 100
    M=int(len(Ek)*P/100)
    figure, axes = subplots(2,2,figsize=(16,8))
## Plot Kinetic Energy Modes ##
    for i in range(size(Ekm,1)):
        if max(Ekm[:,i]>LEGEND_THRESHOLD):
            #axes[0,0].semilogy(t[-M::10],Ekm[-M::10,i],color=colors[i],label=f'$E_{{{i:d}}}$')
            axes[0,0].plot(t[-M::10],Ekm[-M::10,i],color=colors[i],label=f'$E_{{{i:d}}}$')
        else:
            #axes[0,0].semilogy(t[-M::10],Ekm[-M::10,i],color=colors[i])
            axes[0,0].plot(t[-M::10],Ekm[-M::10,i],color=colors[i])

    axes[0,0].grid()
    axes[0,0].set_xlabel(rf'$\Omega_0 t$',fontsize=labelsize,labelpad=labelpadx)
    axes[0,0].set_ylabel(rf'$E_k$',rotation=0,fontsize=labelsize,labelpad=labelpady)
    axes[0,0].set_xlim(t[-M],t[-1])
    axes[0,0].legend(loc='upper right', bbox_to_anchor=(-.25, 1), shadow=True, ncol=1,fontsize=legendfontsize)

## Plot Kinetic Energy Modes above threshold   ##
## Plot Total Kinetic Energy and Knife Edge BC ##
    axes[0,1].semilogy(t[-M::10],Ek[-M::10],'k--',linewidth=2.5,label=f'$E_k$')
    Emax = amax(Ek)
    Emin = amin(Ek)
    for i in range(size(Ekm,1)):
        if max(Ekm[:,i]>PLOT_THRESHOLD):
            axes[0,1].semilogy(t[-M::10],Ekm[-M::10,i],color=colors[i])#,label=f'$E_{{{i:d}}}$')
            if amax(Ekm[-M:,i],axis=0) > Emax:
                Emax = amax(Ekm[-M:,i],axis=0)
            if amin(Ekm[-M:,i],axis=0) < Emin:
                Emin = amin(Ekm[-M:,i],axis=0)
    axes[0,1].grid()
    axes[0,1].set_xlabel(rf'$\Omega_0 t$',fontsize=labelsize,labelpad=labelpadx)
    axes[0,1].set_ylabel(rf'$E_k$',rotation=0,fontsize=labelsize,labelpad=labelpady)
    axes[0,1].set_xlim(t[-M],t[-1])
    axes[0,1].set_ylim(0.1*Emin,100*Emax)

    ax2 = axes[0,1].twinx()
    ax2.semilogy(t[-M::10],w[-M::10],'r-.',label=f'$\omega$')
    ax2.set_ylabel(rf'$\omega$',rotation=0,fontsize=labelsize,labelpad=labelpady)
    linesA, labelsA = axes[0,1].get_legend_handles_labels()
    linesB, labelsB = ax2.get_legend_handles_labels()
    ax2.legend(linesA+linesB,labelsA+labelsB,loc='upper right', bbox_to_anchor=(1.30, 1),
               shadow=True, ncol=1,fontsize=legendfontsize)

## Calculate wMEk using FFT ##
    P  = 80
    M  = int(len(Ek)*P/100)
    dt = t[1]-t[0]
    T  = M*dt     # Period ??
    w0 = 2*pi/T   # Natural Frequency??
    AEk = Ek[-M:].std() # Amplitud of Oscillation
    fftEk  = abs(fft(detrend(Ek[-M:])*blk(M))[:M//2]) # FFT with Blackman filter [array]
    wMEk = w0*fftEk.argmax() # Compute dominant frequency

    #wFFT = min([wMEk,wMEg,wMEw,wMur,wMuw,wMuz])
    wFFT = wMEk
    if wFFT != 0:
        TFFT = 2*pi/wFFT

    # Plot Nperiods periods of Kinetic Energy
    Nperiods = 6
    if Ro>0 and wf>0:
        M=1   #TODO
    elif wFFT>0:
        M  = int(Nperiods*2*pi/(wFFT*dt))
        if M > len(t):
            M=int(len(Ek)*P/100)
    else:
        M=int(len(Ek)*P/100)
## Plot Total Knietic Energy, not log plot ##
    axes[1,0].plot(t[-M:],Ek[-M:])
    axes[1,0].grid()
    axes[1,0].set_xlabel(rf'$\Omega_0 t$',fontsize=labelsize,labelpad=labelpadx)
    axes[1,0].set_ylabel(rf'$E_k$',rotation=0,fontsize=labelsize,labelpad=labelpady)
    axes[1,0].set_xlim(t[-M],t[-1])
    axes[1,0].set_ylim(0.9999*min(Ek[-M:]),1.0001*max(Ek[-M:]))

## Plot FFT ##
    wLim = 1
    AnotationSize = 15
    xPosText = 0.25
    yPosText = 0.92
    ticksize = 12
    # Global Kinetic Energy FFT
    axes[1,1].semilogy(w0*np.arange(len(fftEk)),fftEk,'k-')
    axes[1,1].annotate('$\omega^*$ = {:f}'.format(wMEk), xy=(wMEk, fftEk.max()),
            xycoords='data', xytext=(xPosText,yPosText), textcoords='axes fraction',
            size=AnotationSize, arrowprops=dict(arrowstyle="->"))
    axes[1,1].set_xlabel('$\omega$',fontsize=labelsize,labelpad=labelpadx)
    axes[1,1].set_ylabel('$|\hat{E}_k|$',rotation=0,fontsize=labelsize,labelpad=labelpady)
    axes[1,1].set_xlim(0,wLim)
    axes[1,1].tick_params(labelsize=ticksize)
    if wFFT == 0:
        axes[1,1].set_title(rf'$T = \infty$',fontsize=titlesize)
    else:
        axes[1,1].set_title(rf'$T = {TFFT:f}$',fontsize=titlesize)

    figure.tight_layout()
    savefig(f'{fig_dir:s}Ek{bn:s}',bbox_inches='tight')
    print(f'plotted: {fig_dir:s}Ek{bn:s}')

    print('  Bo  = {:f}'.format(Bo))
    print('  Re  = {:f}'.format(Re))
    print('  Ro  = {:f}'.format(Ro))
    print('  wf  = {:f}'.format(wf))
    print('Gamma = {:f}'.format(gamma))
    print('  T   = {:f}'.format(2*pi/wFFT))
    print('<E_0> = {:f}'.format(average(Ekm[-M:,0])))
    print('std E_0   = {:f}'.format(std(Ekm[-M:,0])))
    print('Delta E_0 = {:f}'.format(max(Ekm[-M:,0])-min(Ekm[-M:,0])))
    print(' ')
    return None

def main(f):
    FIG_DIR = get_figdir(f)
    os.makedirs(FIG_DIR,exist_ok=True)
    fields = get_fields(f)
    plotter(fields,FIG_DIR)
    return None

if __name__ == '__main__':
    G = natsort.realsorted(glob.glob(sys.argv[1]))
    # NOTE G = chec_if_plotted()?
    p = mp.Pool(processes=NPROCS)
    p.map(main,G)
