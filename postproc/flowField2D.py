#!/usr/bin/env python
import sys,os
import argparse
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial
from glob import glob
from myPlots import * # It also imports pylab
from myReader import * # It also imports numpy

def get_args():
  class myFormatter(argparse.RawDescriptionHelpFormatter,
          argparse.ArgumentDefaultsHelpFormatter):
          pass
  parser = argparse.ArgumentParser(description=main.__doc__,
          formatter_class=myFormatter)
  parser.add_argument('DataFilePath',
          help='Path to input file, restart or vtk')
  parser.add_argument('Problem',
          help='Name of the problem to plot: kedgeTop/freeSurfTop')
  parser.add_argument('Method',
          help='Method used for DNS: FD/Spectral')
  parser.add_argument('Fields',
          help='List of fields to plot')
  parser.add_argument('Clips',type=list,
          help='List of clipping values for each field')
  parser.add_argument('-s','--size',dest='figBaseSize', type=int,
          help='Base size of the figure',default=4)
  parser.add_argument('-g','--gamma-disable',dest='gammaOpt',action='store_false',
          help='Aspect ratio of figure: Height/Lenth',default=True)
  parser.add_argument('-m','--max',dest='gma', type=float,
          help='Global absolute maximum. gma = max(abs(f)) where f is the \
          field')
  #parser.add_argument('-c','--create',action='store_true',
          #help='Aspect ratio of figure: Height/Lenth',default=False)
  args = parser.parse_args()

  fields = args.Fields.split(',')
  clips = ''.join(args.Clips).split(',')
  clips = list(map(float, clips))

  return args.DataFilePath,args.Problem,args.Method,fields,clips,\
      args.figBaseSize,args.gammaOpt,args.gma

def get_files_to_plot(DataFilePath):
  drecs = []
  if '*' in DataFilePath:
    drecs = glob(DataFilePath) # list of input restarts, path included
    if not drecs:
      print('\nNO FILES FOUND MATCHING THE PATTERN:  ', DataFilePath)
      print('--> QUITTING')
      sys.exit(1)
  else:
    drecs.append(DataFilePath)
    if not os.path.exists(DataFilePath):
      print('\nFILE ',DataFilePath,' NOT FOUND')
      print('--> QUITTING')
      sys.exit(1)
  drecs.sort()
  return drecs

def mycf(X,Y,field,clip,fn=1,fb=4,fgamma=1.0,Nell=33,Nred=3,gma=None):
  if not gma: # Calculate abolute maximum if gma not specified
    gma = abs(field).max()
  ell = clip*pylab.linspace(-1,1,Nell)*gma
  red = pylab.linspace(clip*gma,gma,Nred)
  blu = -red[::-1]
  f,a = no_ax_fig(fn,fb,fgamma)
  mycontourf(X,Y,field,lc='#777777',lw=0.1,levels=blu,cmap=myBlues)
  mycontourf(X,Y,field,lc='#777777',lw=0.1,levels=ell,cmap=mycm15)
  mycontourf(X,Y,field,lc='#777777',lw=0.1,levels=red,cmap=myReds)
  #if field.min() < 0 and field.max() > 0:
  #    contour(X,Y,field,levels=[0],linestyles='-',colors='#777777')
  return None

def main(f,problem,method,fields,clips,fb,gammaOpt,gma):
    '''
    Generates a figure of size (figBaseSize,Gamma*figBaseSize) with 1 axes
    The axes and the box around the figure are not visible
    '''

    outDir = 'fig/'
    try:
        os.mkdir(outDir)
    except:
        pass

    bn = os.path.basename(f) # basename
    ext = os.path.splitext(f)[-1] # extension

    if method.lower() == 'fd':
      if ext == '.vtk':
        print('\nFD codes do not work with vtk files')
        sys.exit(1)
      elif problem.lower() == 'kedgetop':
        R,Z,d = read_kedgeTop_restart(f)
      elif problem.lower() == 'freesurftop':
        R,Z,d = read_freeSurfTop_restart(f)
    elif method.lower() == 'spectral':
      if ext is not '.vtk':
        print('\nSpectral codes not work with vtk files')
        sys.exit(1)
      elif problem.lower() == 'kedgetop':
        R,Z,d = read_vtk(f)
      elif problem.lower() == 'freesurftop':
        print('\nNot coded yet!')
        sys.exit(1)

    for i in range(len(fields)):
      if gammaOpt:
        # Obtain parameters from basename
        tokens = find_tokens(bn)
        params = {token:float(parse_token(bn,token)) for token in tokens}
        mycf(R,Z,d[fields[i]],clips[i],fb=fb,fgamma=params['Gamma'],gma=gma)
      else:
        mycf(R,Z,d[fields[i]],clips[i],fb=fb,gma=gma)
      plt.savefig(outDir+fields[i]+'_'+bn.split('.')[0]+'.png')
      plt.close()
    return None

def check_input_args(probName, method, fields, clips):
  allowedProblems = ['kedgetop','freesurftop']
  allowedMethods  = ['fd','spectral']

  if probName.lower() not in allowedProblems:
    print('\nThe problem name is not correct.')
    print('\nAvailable names: ', *allowedProblems,sep="\n -- ")
    sys.exit(1)
  elif method.lower() not in allowedMethods:
    print('\nThe method is not correct.')
    print('\nAvailable methods: ', *allowedMethods,sep="\n -- ")
    sys.exit(1)
  elif len(fields) is not len(clips):
    print('\nThe lists Fields and Clips must have the same lenth.')
    print('\nOne clipping value per field!')
    sys.exit(1)
  elif max(clips)>=1 or min(clips)<=0:
    print('\nClips out of range.')
    print('\nThe Clips must be in (0,1).')
    sys.exit(1)
  return None

if __name__ == '__main__':
  DataFilePath, probName, method, fields, clips, \
   figBaseSize, gammaOpt, gma = get_args()

#  if create: os.system("./myParaview < pv_input &")

  check_input_args(probName, method, fields, clips)

  drecs = get_files_to_plot(DataFilePath)

  NPROCS = 4
  if NPROCS > 1:
    mkl_set_num_threads(NPROCS)
    pool = mp.Pool(processes=NPROCS)
    D  = partial(main,problem=probName,method=method,fields=fields,
                      clips=clips,fb=figBaseSize,gammaOpt=gammaOpt,gma=gma)
    DD = pool.map(D,drecs)
  else:
    for drec in drecs:
      main(drec)

