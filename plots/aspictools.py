#!/usr/bin/python
import numpy as np
import aspicio as aio


def log_energy_reheat_ingev(lnRhoReh):
    lnMpinGeV= 42.334
    logErehGeV = 0.25*(lnRhoReh + 4*lnMpinGeV)/np.log(10)
    return logErehGeV

def r_epstwo_zero(ns):
    r = 8*(1-ns)
    return r

def ns_epstwo_zero(r):
    ns = 1-r/8.
    return ns

def r_epstwo_twoepsone(ns):
    r = 4*(1-ns)
    return r

def ns_epstwo_twoepsone(r):
    ns = 1-r/4.
    return ns

def separatrix_rns(rmin,rmax,npoints=50):
    r = np.linspace(rmin,rmax,npoints)
    ns1 = ns_epstwo_twoepsone(r)
    ns2 = ns_epstwo_zero(r)    
    return r,ns1,ns2


def separatrix_nsr(nsmin,nsmax,npoints=50):
    ns = np.linspace(nsmin,nsmax,npoints)
    r1 = r_epstwo_twoepsone(ns)
    r2 = r_epstwo_zero(ns)
    return ns,r1,r2


def separatrix_epstwo_twoepsone(xmin,xmax,npoints=50):
    eps2 = np.exp(np.linspace(np.log(xmin),np.log(xmax),npoints))
    eps1 = 0.5*eps2
    return eps2,eps1



def contour_plot(ax,file1,col1,file2,col2,alpha,depth=None,style=None,log=None,**kwargs):

    xy1s = list(aio.read_contour(file1))
    xy2s = list(aio.read_contour(file2))


    if log:
        ax.fill(xy1s[0],10**xy1s[1],col1,xy2s[0],10**xy2s[1],col2
                ,zorder=depth,alpha=alpha,linestyle=style,**kwargs)
    else:
        ax.fill(xy1s[0],xy1s[1],col1,xy2s[0],xy2s[1],col2
                ,zorder=depth,alpha=alpha,linestyle=style,**kwargs)

    return


def model_prediction(ax,namei,fooi,cmap,marker,alpha
                     ,ilab=None,txtpos=None,cstyle=None):
    depth = 6
    scmodel = ax.scatter(fooi[0],fooi[1],c=fooi[3],zorder=depth,marker=marker
                         ,cmap=cmap,alpha=alpha)

    if ilab is not None:
        ax.annotate(namei, xy=(fooi[0][ilab],fooi[1][ilab]),  xycoords='data',
                    xytext=txtpos, textcoords='offset points',
                    zorder=depth, color='black',
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle=cstyle))
        print 'Ereh = ',fooi[3][ilab]
    return scmodel


def model_prediction_efold(ax,namei,fooi,cmap,marker,alpha
                     ,ilab=None,txtpos=None,cstyle=None):
    depth = 6
    scmodel = ax.scatter(fooi[0],fooi[1],c=fooi[2],zorder=depth,marker=marker
                         ,cmap=cmap,alpha=alpha)

    if ilab is not None:
        ax.annotate(namei, xy=(fooi[0][ilab],fooi[1][ilab]),  xycoords='data',
                    xytext=txtpos, textcoords='offset points',
                    zorder=depth, color='black',
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle=cstyle))
        print 'DeltaN* = ',fooi[2][ilab]
    return scmodel


