#!/usr/bin/python
import numpy as np
import aspicio as aio
import aspicfigs as asf
import matplotlib.colors as mplc

def reduced_cmap(numcolors=32, name='reducecolor',
               colors=['black','darkred','darkorange','orange','yellow','white']):

    
    cmap = mplc.LinearSegmentedColormap.from_list(name=name, 
                                             colors = colors,N=numcolors)
    return cmap


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


def r_epstwo_is_alpha_epsone(ns,alpha):
    r = 16*(1-ns)/(2+alpha)
    return r

def ns_epstwo_is_alpha_epsone(r,alpha):
    ns = 1-(2+alpha)*r/16
    return r



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



def separatrix_rns_alpha(alpha,rmin,rmax,npoints=50):
    r = np.linspace(rmin,rmax,npoints)
    ns1 = ns_epstwo_is_alpha_epsone(r,alpha)
    return r,ns1


def separatrix_nsr_alpha(alpha,nsmin,nsmax,npoints=50):
    ns = np.linspace(nsmin,nsmax,npoints)
    r1 = r_epstwo_is_alpha_epsone(ns,alpha)
    return ns,r1


def separatrix_epstwo_twoepsone(xmin,xmax,npoints=50):
    return separatrix_epstwo_is_alpha_epsone(alpha=2.0,xmin=xmin,xmax=xmax,npoints=npoints)


def separatrix_epstwo_is_alpha_epsone(alpha,xmin,xmax,npoints=50):
    eps2 = np.exp(np.linspace(np.log(xmin),np.log(xmax),npoints))
    eps1 = eps2/alpha
    return eps2,eps1



def contour_plot(ax,file1,col1,file2,col2,alpha,depth=None,style=None,log=None,**kwargs):

    xy1s = list(aio.read_contour(file1))
    xy2s = list(aio.read_contour(file2))


    if log:
        ax.fill(xy1s[0],10**xy1s[1],col1,xy2s[0],10**xy2s[1],col2
                ,zorder=depth,alpha=alpha,linestyle=style,edgecolor='black',**kwargs)
    else:
        ax.fill(xy1s[0],xy1s[1],col1,xy2s[0],xy2s[1],col2
                ,zorder=depth,alpha=alpha,linestyle=style,edgecolor='black',**kwargs)

    return


def model_prediction(ax,namei,fooi,cmap,marker,alpha,
                     ilab=None,txtpos=None,cstyle=None,colindex=3,linewidth=0.5,markersize=16,
                     edgecolors='black'):
    depth = 6
    scmodel = ax.scatter(fooi[0],fooi[1],s=markersize,c=fooi[colindex],zorder=depth,marker=marker
                         ,cmap=cmap,alpha=alpha,linewidths=linewidth,edgecolors=edgecolors)

    
    if ilab is not None:
        ax.annotate(namei, xy=(fooi[0][ilab],fooi[1][ilab]),  xycoords='data',
                    xytext=txtpos, textcoords='offset points',
                    zorder=depth, color='black',
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle=cstyle))
        print('Ereh = ',fooi[3][ilab])

   
    return scmodel


def model_prediction_efold(ax,namei,fooi,cmap,marker,alpha
                     ,ilab=None,txtpos=None,cstyle=None):
    depth = 6
    scmodel = ax.scatter(fooi[0],fooi[1],c=fooi[2],zorder=depth,marker=marker
                         ,cmap=cmap,alpha=alpha,edgecolors='black')

    if ilab is not None:
        ax.annotate(namei, xy=(fooi[0][ilab],fooi[1][ilab]),  xycoords='data',
                    xytext=txtpos, textcoords='offset points',
                    zorder=depth, color='black',
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle=cstyle))
        print('DeltaN* = ',fooi[2][ilab])
    return scmodel


def model_prediction_filtered(fig,ax,xybounds,xyzoomrange,namei,fooi,scmap,marker,scalpha,
                              semilog=True,threshlabel=None,threshscat=None,
                              ilab=None,txtpos=None,cstyle=None):

    depth = 6

    if threshlabel is None:
        threshlabel = 0.45

    if threshscat is None:
        threshscat = 0.1
        
    xycscat, sannotate, xyannotate, xylabel, labangles, countplot, countrange = asf.filter_model_data(
        fig.get_size_inches(),xybounds,xyzoomrange,fooi,semilog,threshscat,threshlabel)
    
    scmodel = asf.add_model_scatter(ax,xycscat,scmap,marker=marker,alpha=scalpha,depth=depth)


    
    if ilab is not None:
        par = fooi[ilab]
        ax.annotate(namei, xy=(par[0].getvalue(),par[1].getvalue()),  xycoords='data',
                    xytext=txtpos, textcoords='offset points',
                    zorder=depth, color='black',
                    arrowprops=dict(arrowstyle="->",
                                    connectionstyle=cstyle))
        print('Ereh = ',par[3].getvalue())

    return scmodel
