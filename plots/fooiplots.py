#!/usr/bin/python
import numpy as np
import re as re
import argparse
import glob
import aspicio as aio
import aspicfigs as asf

def swap_in(alist,f,t):
    alist[f], alist[t] = alist[t], alist[f]


def tryint(s):
    try:
        return int(s)
    except:
        return s
    
def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [ tryint(c) for c in re.split('([0-9]+)', s) ]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)




parser = argparse.ArgumentParser()

parser.add_argument("dir", help="data file directory")
parser.add_argument("name", help="model name")
parser.add_argument("type", help="'powerlaw' or 'slowroll' data file")
parser.add_argument("--ste", action="store_true", help="draw STE separatrix")
parser.add_argument("--fillste", action="store_true", help="color fill the STE regions")
parser.add_argument("--contour", action="store_true",help="draw one- and two-sigma contours")
parser.add_argument("--fillcontour", action="store_true",help="color fill one- and two-sigma contours")
parser.add_argument("--semilog", action="store_true",help="vertical axis in logarithmic scale")
parser.add_argument("--wreh", type=str, help="value of wreh if non zero, in laTeX format, e.g. '$-1/3$'")

parser.add_argument("--decimalfmt", type=str, help="format to use to display decimal parameter values, e.g. '%%.2g'")
parser.add_argument("--threshscat", type=float, help="minimal separation between points (in inches)")
parser.add_argument("--modulolabel", type=int, help="skip labels using this modulo conditional")
parser.add_argument("--tiltlabel", type=float, help="add an additional tilt angle to labels, in degrees")
parser.add_argument("--threshlabel", type=float, help="minimal separation between labels (in inches)")
parser.add_argument("--movielabel", type=float, help="add a tilt angle beetween two consecutive labels, in degrees")
parser.add_argument("--arrowlength", type=float, help="set the size of the labels'arrow, in pts")
parser.add_argument("--xyinibounds", type=float, nargs="+", help="initial values of the window boundary: xmin xmax ymin ymax")
parser.add_argument("--nozoomout", action="store_true", help="forbids zooming out")
parser.add_argument("--zoomplot", type=int, help="zoom out till the number of points visible in the plot \
reached the input value, or till the number of steps reaches --zoomsteps or till the window reaches --xyhardbounds")
parser.add_argument("--zoomsteps", type=int, help="number of zooming out steps allowed")
parser.add_argument("--xyhardbounds", type=float, nargs="+", help="absolute minimal and maximal \
values of all axis: min(xmin) max(xmin) min(xmax) max(xmax) min(ymin) max(ymin) min(ymax) max(ymax)")
parser.add_argument("--zoomrange", type=int, help="zoom out till the number of points visible \
in --xyzoomrange range reaches the input value")
parser.add_argument("--xyzoomrange", type=float, nargs="+", help="count the number of model's \
points only within the sub-window: xmin xmax ymin ymax")

pargs = parser.parse_args()


aio.read_modeldict('modeldict.txt')
aio.read_paramdict('paramdict.txt')
model = aio.model(name=pargs.name)

if pargs.type == 'powerlaw':

    swapxy = False

    if pargs.xyinibounds is not None:
        xyinibounds = tuple(pargs.xyinibounds)
    else:
        xyinibounds = (0.93,1.005,0.0005,0.2)

    if pargs.xyhardbounds is not None:
        xyhardbounds  = [(pargs.xyhardbounds[0],pargs.xyhardbounds[1]),
                         (pargs.xyhardbounds[2],pargs.xyhardbounds[3]),
                         (pargs.xyhardbounds[4],pargs.xyhardbounds[5]),
                         (pargs.xyhardbounds[6],pargs.xyhardbounds[7])]
    else:
        xyhardbounds  = [(0.88,0.945),(0.98,1.2),(1e-30,1e-3),(2e-3,0.8)]

    if pargs.xyzoomrange is not None:
        xyzoomrange = tuple(pargs.xyzoomrange)
    else:
        xyzoomrange = (0.95,0.97,1e-20,0.8)

    
    contourfiles=('data/contour_ns_logr_level_0.dat'
                 ,'data/contour_ns_logr_level_1.dat')

    figprefix = 'nsr_'+ model.getname()

elif pargs.type == 'slowroll':

    swapxy = True

    if pargs.xyinibounds is not None:
        xyinibounds = tuple(pargs.xyinibounds)
    else:
        xyinibounds = (0.005,0.08,9e-6,1e-2)

    if pargs.xyhardbounds is not None:
        xyhardbounds  = [(pargs.xyhardbounds[0],pargs.xyhardbounds[1]),
                         (pargs.xyhardbounds[2],pargs.xyhardbounds[3]),
                         (pargs.xyhardbounds[4],pargs.xyhardbounds[5]),
                         (pargs.xyhardbounds[6],pargs.xyhardbounds[7])]
    else:
        xyhardbounds  = [(-0.1,0.02),(0.05,0.15),(1e-40,1e-3),(2e-3,0.8)]

    if pargs.xyzoomrange is not None:
        xyzoomrange = tuple(pargs.xyzoomrange)
    else:
        xyzoomrange = (0.03,0.05,1e-20,0.5)

    
    contourfiles=('data/contour_sr2_sr1_level_0.dat'
                 ,'data/contour_sr2_sr1_level_1.dat')

    figprefix = 'eps_'+ model.getname()

else:

    raise Exception("type unknown: only 'powerlaw' and 'slowroll' are currently recknownized!")



    

if pargs.threshscat is not None:
    myThreshScat = pargs.threshscat
else:
    myThreshScat = 0.1

    

if pargs.threshlabel is not None:
    myThreshLabel = pargs.threshlabel
else:
    myThreshLabel = 0.8


if pargs.modulolabel is not None:
    myModLabel = pargs.modulolabel
else:
    myModLabel = 1


if pargs.tiltlabel is not None:
    myTiltLabel = pargs.tiltlabel*np.pi/180.0
else:
    myTiltLabel = 0.0


if pargs.movielabel is not None:
    myMovieLabel = pargs.movielabel*np.pi/180.0
else:
    myMovieLabel = None


    
fileprefix = pargs.dir+'/'+model.getname()+'_'+pargs.type+'_'

aspicfiles = glob.glob(fileprefix+'*.dat')
sort_nicely(aspicfiles)

for i in range(len(aspicfiles)):
    datafile = aspicfiles[i]
    figfile = figprefix + '_' + str(i) + '.eps'

    print('Reading data',datafile)
    fixedParams, varParams = aio.read_aspicfile(datafile)

    if swapxy:
        for par in varParams:
            swap_in(par,0,1)
    
    print('Generating figure',figfile)
    asf.create_figure(model,fixedParams,varParams,contourfiles,figfile,xyinibounds,
                      xyhardbounds,xyzoomrange,pargs,threshscat=myThreshScat,
                      threshlabel=myThreshLabel, modlabel=myModLabel, tiltlabel=myTiltLabel,
                      movielabel=myMovieLabel)

    print

