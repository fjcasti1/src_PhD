#!/usr/bin/env python
import sys,os
import argparse
import matplotlib.pyplot as plt
import multiprocessing as mp
from functools import partial
from glob import glob
from myPlots import * # It also imports pylab

def get_args():
    class myFormatter(argparse.RawDescriptionHelpFormatter,
            argparse.ArgumentDefaultsHelpFormatter):
            pass
    parser = argparse.ArgumentParser(description=main.__doc__,
            formatter_class=myFormatter)
    parser.add_argument("RestartFilesPath",
            help="Path to restart files. If there is a glob, it requires to be \
            between single quotes. VTK files automatically ignored")
    parser.add_argument("Fields",
            help="List of fields to generate")
    parser.add_argument("Dimensions",
            help="List of dimensions Nr,Nz,Nt of the physical uniform grid. \
            Nt = 0 for Axisymmetric cases")
    parser.add_argument("-m","--mode", type=int,
            help="[ +m: 0+m | -m: m | +1000: all modes | -1000: all modes but 0 ]",
            default=1000)
    args = parser.parse_args()


    fields = args.Fields.split(',')
    if fields[0]=='all':
        fields = ['vr','vt','vz',\
                  'wr','wt','wz',
                  'rv','ke','he',
                  'sf','vx','vy']

    dims = ''.join(args.Dimensions).split(',')
    dims = list(map(int, dims))
    return args.RestartFilesPath,fields,dims,args.mode

def get_files_to_use(DataFilePath):
    drecs = []
    if '*' in DataFilePath:
        # list of input restarts, path included
        drecs = [drec for drec in glob(DataFilePath) if '.vtk' not in drec]
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

def mycf(X,Y,field,clip,fn=1,fb=4,fgamma=1,Nell=33,Nred=3):
    gma = abs(field).max()
    ell = clip*pylab.linspace(-1,1,Nell)*gma
    red = pylab.linspace(clip*gma,gma,Nred)
    blu = -red[::-1]
    f,a = no_ax_fig(fn,fb,fgamma)
    mycontourf(X,Y,field,lc='#777777',lw=1,levels=blu,cmap=myBlues)
    mycontourf(X,Y,field,lc='#777777',lw=0,levels=ell,cmap=mycm15)
    mycontourf(X,Y,field,lc='#777777',lw=1,levels=red,cmap=myReds)
    #if field.min() < 0 and field.max() > 0:
    #    contour(X,Y,field,levels=[0],linestyles='-',colors='#777777')
    return None

def create_input(drecs,fields,dims,mode):
    # Field dictionary options
    fd = {"vr":0,"vt":0,"vz":0,
          "wr":0,"wt":0,"wz":0,
          "rv":0,"ke":0,"he":0,
          "sf":0,"vx":0,"vy":0}
    for field in fields:
        fd[field] = 1

    pv_input = f"""{dims[0]:d} ! nrp
{dims[1]:d} ! nzp
{dims[2]:d} ! ntp
{mode:d}    ! i0mode [ +m: 0+m | -m: m | +1000: all modes | -1000: all modes but 0 ]
   {fd['vr']:d}  {fd['vt']:d}  {fd['vz']:d} \
 {fd['wr']:d}  {fd['wt']:d}  {fd['wz']:d} \
 {fd['rv']:d}  {fd['ke']:d}  {fd['he']:d} \
 {fd['sf']:d}  {fd['vx']:d}  {fd['vy']:d}  [0 not plotted, 1 plotted]
! vr vt vz wr wt wz rv ke he sf vx vy
! ------------------------------------------------------------------------------
! WARNING: be sure the subdirectory movie/ exists!!!
! ------------------------------------------------------------------------------
! list of files to be processed:          !a blank line interrupts the execution
"""
    for drec in drecs:
        pv_input+=drec+'\n'
    return pv_input

def main(f,fields,clips,fb,Gamma):
    '''
    Generates an input file for the Fortran vtk generator and runs it
    generating the vtk files.
    '''

    os.system("./bin/myParaview < pv_input &")


    return None

if __name__ == '__main__':
    RestartFilesPath,fields,dims,mode = get_args()

    drecs = get_files_to_use(RestartFilesPath)

    pv_input = create_input(drecs,fields,dims,mode)

    with open('pv_input','w') as f:
        f.write(pv_input)
        f.close()
    os.system("./bin/myParaview < pv_input")


#    NPROCS = 4
#    if NPROCS > 1:
#      mkl_set_num_threads(NPROCS)
#      pool = mp.Pool(processes=NPROCS)
#      D  = partial(main,fields=fields,clips=clips,
#              fb=figBaseSize,Gamma=Gamma)
#      DD = pool.map(D,drecs)
#    else:
#      for drec in drecs:
#        main(drec)

