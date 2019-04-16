#!/usr/bin/python
import numpy as np
import aspicio as aio
import aspictools as ast
import matplotlib.colors as mplc


def aspic_cmap(numcolors=32, name='aspicolor',
               colors=['black','darkred','darkorange','orange','yellow','white']):

    
    cmap = mplc.LinearSegmentedColormap.from_list(name=name, 
                                             colors = colors,N=numcolors)
    return cmap


def sci_notation(num, decimal_digits=1, precision=None, exponent=None, expmin=3, decfmt='%.2g'):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """

    if num == 0.0:
        return r'$%s$' % float(decfmt % num)

    if not exponent:
        exponent = int(np.floor(np.log10(abs(num))))        
    coeff = np.round(num / float(10**exponent), decimal_digits)
    if not precision:
        precision = decimal_digits

    if np.abs(exponent) > expmin:
        return r"${0:.{2}f}\times10^{{{1:d}}}$".format(coeff, exponent, precision)
    else:
        return r'$%s$' % float(decfmt % num)


def get_axsize(ax,fig):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height



def get_xydatabounds(allpars,enlargeby=0.1):

    bounds =(+np.inf,-np.inf,+np.inf,-np.inf)
    
    for par in allpars:

        bounds = (min(bounds[0],par[0].getvalue()),
                  max(bounds[1],par[0].getvalue()),
                  min(bounds[2],par[1].getvalue()),
                  max(bounds[3],par[1].getvalue()))


    if enlargeby != 0:
        bounds = (bounds[0] - enlargeby*np.abs(bounds[0]),
                  bounds[1] + enlargeby*np.abs(bounds[1]),
                  bounds[2] - enlargeby*np.abs(bounds[2]),
                  bounds[3] + enlargeby*np.abs(bounds[3]))

    
    return bounds
        


def get_zoomdirections(xycurrent,xytarget,semilog = False):

    angles = [0,1]
    dists = [0,1]

    if semilog:
        xycurrent = [xycurrent[0],xycurrent[1],np.log10(xycurrent[2]),np.log10(xycurrent[3])]
        xytarget = [xytarget[0],xytarget[1],np.log10(xytarget[2]),np.log10(xytarget[3])]
    
    angles[0] = np.arctan2(xytarget[2]-xycurrent[2], xytarget[0]-xycurrent[0])
    angles[1] = np.arctan2(xytarget[3]-xycurrent[3], xytarget[1]-xycurrent[1])

    dists[0] = np.sqrt((xytarget[2]-xycurrent[2])**2 + (xytarget[0]-xycurrent[0])**2)
    dists[1] = np.sqrt((xytarget[3]-xycurrent[3])**2 + (xytarget[1]-xycurrent[1])**2)

    return angles, dists



def get_zoomxybounds(xybounds,angles,steps,semilog = False):

    xyzoom = [0,1,2,3]

    xyzoom[0] = xybounds[0] + steps[0] * np.cos(angles[0])
    xyzoom[1] = xybounds[1] + steps[1] * np.cos(angles[1])
    if semilog:
        xyzoom[2] = 10**(np.log10(xybounds[2]) + steps[0] * np.sin(angles[0]))
        xyzoom[3] = 10**(np.log10(xybounds[3]) + steps[1] * np.sin(angles[1]))
    else:
        xyzoom[2] = xybounds[2] + steps[0] * np.sin(angles[0])
        xyzoom[3] = xybounds[3] + steps[1] * np.sin(angles[1])

    return xyzoom


def check_hardxybounds(xy,allminmax):

    xycheck = [0,1,2,3]

    for i in range(len(allminmax)):
        minmax = allminmax[i]
        if xy[i] < minmax[0]:
            xycheck[i] = minmax[0]
        elif xy[i] > minmax[1]:
            xycheck[i] = minmax[1]
        else:
            xycheck[i] = xy[i]


    return xycheck


def filter_model_data(inches,xybounds,xyrange,allpars,semilog=False,threshscat=0.1,
                      threshlabel=1.0,modlabel=1,decfmt='%.2g'):

    xmin = xybounds[0]
    xmax = xybounds[1]
    ymin = xybounds[2]
    ymax = xybounds[3]    
    
    xscat = []
    yscat = []
    cscat = []

    xyannotate = []
    sannotate = []
    labangles = []
    

    parval = None
    xysavscat = (None,None)
    xysavlab = (None,None)
    labangle = 0.0
    countlabel = 0
    
    count = 0
    
    plotscat=True

    countrange = 0
    inrange = False

    for par in allpars:
        
        xinch = (par[0].getvalue()-xmin) * inches[0]/(xmax-xmin)
        if not semilog:
            yinch = (par[1].getvalue()-ymin) * inches[1]/(ymax-ymin)
        else:
            yinch = (np.log10(par[1].getvalue())-np.log10(ymin))*inches[1]/(np.log10(ymax)-np.log10(ymin))
            
        if xinch < 0 or yinch < 0 or xinch > inches[0] or yinch > inches[1]:
            continue
        
        if len(par) < 5:
            plotscat = True


        elif par[4].getvalue() <> parval:
            if count <> 0:
                labangles.append(labangle)
            parval = par[4].getvalue()
            count = 0
            labangle = 0.0
            if not_overlap(xysavlab,(xinch,yinch),threshlabel):
                xysavlab = (xinch,yinch)
                xyannotate.append((par[0].getvalue(),par[1].getvalue()))

                if countlabel % modlabel == 0:
                    sannotate.append(par[4].getlabel() + r'$=$'+sci_notation(par[4].getvalue()
                                                                             ,decfmt=decfmt))
                else:
                    sannotate.append(None)
                    
                countlabel += 1
                plotscat = True


                
            else:
                plotscat = False

#was elif (!?)                
        if par[4].getvalue() == parval and plotscat:

#unnecessary?            if xysavlab <> (None,None):
            count += 1
            labangle = (labangle*(count-1)
                        + np.arctan2(yinch - xysavlab[1],xinch - xysavlab[0]))/count               
                

        if plotscat and not_overlap(xysavscat,(xinch,yinch),threshscat):
            xysavscat = (xinch,yinch)
            xscat.append(par[0].getvalue())
            yscat.append(par[1].getvalue())
            cscat.append(float(decfmt % ast.log_energy_reheat_ingev(par[3].getvalue())))

            
            inrange = par[0].getvalue() > xyrange[0] and par[0].getvalue() < xyrange[1] \
                      and par[1].getvalue() > xyrange[2] and par[1].getvalue() < xyrange[3]
            
            if inrange:
                countrange += 1


#    print(len(xyannotate), len(labangles), len(sannotate))
#    print("test",xyannotate)
#    print("label",labangles)
                
#by construction the for loop may miss appending the very last label            
    if len(xyannotate) == len(labangles)+1:
        labangles.append(labangle)



        
    return (xscat,yscat,cscat), sannotate, xyannotate, \
           (par[0].getlabel(),par[1].getlabel()), \
           labangles,len(cscat),countrange




def not_overlap(xy1,xy2,thresh):
    if any(x is None for x in xy1) or any(x is None for x in xy2):
        return True
   
    dist = np.sqrt((xy2[0]-xy1[0])**2 + (xy2[1]-xy1[1])**2)
    if dist==0 or dist <= thresh:
        return False
    else:
        return True


    
def add_model_scatter(ax,xycdata,cmap,marker,alpha,depth=1):
    scmodel = ax.scatter(xycdata[0],xycdata[1],c=xycdata[2],zorder=depth,marker=marker
                         ,cmap=cmap,alpha=alpha,clip_on=True)
    return scmodel



def add_model_annotate(ax,s,xy,xytextpts=None,cstyle=None,depth=3):

    if len(s) <> len(xy):
        raise Exception("string and (x,y) of unequal length!")

    mingap = 10
    
    if xytextpts is None:
        xtext = mingap
        ytext = mingap

        
    elif len(xytextpts) <> len(xy):
        print('arg',len(xytextpts),len(xy),len(s))
        raise Exception("(x,y) and (xtext,ytext) of unequal length!")
        
    
    for i in range(len(s)):

        if xytextpts is not None:
            xtext,ytext = xytextpts[i]

        if s[i] is not None:
            ax.annotate(s[i],xy[i],xycoords='data',
                        xytext=(xtext,ytext),textcoords='offset points',
                        fontsize=8, zorder=depth, color='black',
                        clip_on=False, annotation_clip=False, horizontalalignment = 'center',
                        bbox=dict(boxstyle="roundtooth", fc="w", color='slategray'),
                        arrowprops=dict(arrowstyle="->",color='slategray',
                                        connectionstyle=cstyle))


def add_contour_plot(ax,file1,col1,file2,col2,alpha,depth=None,style=None,log=None,**kwargs):

    xy1s = list(aio.read_contour(file1))
    xy2s = list(aio.read_contour(file2))


    if log:
        ax.fill(xy1s[0],10**xy1s[1],col1,xy2s[0],10**xy2s[1],col2
                ,zorder=depth,alpha=alpha,linestyle=style,closed=False,**kwargs)
    else:
        ax.fill(xy1s[0],xy1s[1],col1,xy2s[0],xy2s[1],col2
                ,zorder=depth,alpha=alpha,linestyle=style,closed=False,**kwargs)

    return


    
def add_model_name(ax,string):  
    ax.text(0.97,0.97,string
            , bbox=dict(facecolor='white', alpha=1)
            ,horizontalalignment='right',verticalalignment='top',transform=ax.transAxes)



def add_model_fixparams(ax,pars):
    yfix = 0.97
    yshift = 0.04
    for i in range(len(pars)):
#        string = pars[i].getlabel() + r'$ = %s$' % float('%.2g' % pars[i].getvalue())
        string = pars[i].getlabel() + r'$=$'+sci_notation(pars[i].getvalue())
        ax.text(0.02,yfix,string,horizontalalignment='left',verticalalignment='top'
                , transform=ax.transAxes,bbox=dict(boxstyle="square", pad=0.0, fc="w", color='w'))
        yfix = yfix - yshift



def create_figure(model, allfixed, allpars, contourfile, figfile, xyinibounds, xyhardbounds, xyzoomrange,
                  args, zoomplot=10, zoomrange=0, zoomsteps=10, dpi=150, threshscat=0.1,
                  threshlabel=1.0,modlabel=1,tiltlabel=0.0,movielabel=None):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    filetype = args.type
    fillste = args.fillste
    ste = args.ste
    
    fillcontour = args.fillcontour
    contour = args.contour

    semilog = args.semilog
    wrehbar = args.wreh

    if args.decimalfmt is not None:
        decfmt = args.decimalfmt
    else:
        decfmt = '%.2g'

    if args.zoomplot is not None:
        zoomplot = args.zoomplot

    if args.zoomsteps is not None:
        zoomsteps = args.zoomsteps

    if args.zoomrange is not None:
        zoomrange = args.zoomrange


    scmap = aspic_cmap()
    scalpha = 1

    fig, ax = plt.subplots()
    fig.set_dpi(dpi)

    xydatabounds = get_xydatabounds(allpars)
    
    zoomcount = 0
    xybounds = xyinibounds
    zangles, zdists = get_zoomdirections(xybounds,xydatabounds,semilog)

    zsteps = []
    for dist in zdists:
        zsteps.append(dist/zoomsteps)

        
    while zoomcount < zoomsteps + 1:
        xycscat, sannotate, xyannotate, xylabel, labangles, countplot, countrange = filter_model_data(
            fig.get_size_inches(),xybounds,xyzoomrange,allpars,semilog,threshscat,
            threshlabel,modlabel,decfmt)
        
        if args.nozoomout is not None and args.nozoomout:
            break

        
        if countrange < zoomrange:
            xybounds = get_zoomxybounds(xybounds,zangles,zsteps,semilog)
            xybounds = check_hardxybounds(xybounds,xyhardbounds)
        elif countplot < zoomplot:
            xybounds = get_zoomxybounds(xybounds,zangles,zsteps,semilog)
            xybounds = check_hardxybounds(xybounds,xyhardbounds)
        else:
            break
        
        zoomcount += 1

    if countplot == 0 and countrange == 0:
        print "create_figure: cannot find data within zoom limits!"
        return

    ax.set_autoscale_on(False)
    ax.axis(xybounds)

    arc = +0.1


    if semilog:
        ax.semilogy()
        

    xyptoffset = []
    rptoffset = 40
    
    addangle = np.pi + tiltlabel
    
    for i in range(len(labangles)):

        xptoffset = rptoffset*np.cos(labangles[i]+addangle)
        yptoffset = rptoffset*np.sin(labangles[i]+addangle)

        xyptoffset.append((xptoffset,yptoffset))
            
        if movielabel is not None:
            addangle = addangle + movielabel
    

#Swartchz Terro-Escalante regions
    if fillste or ste:

        npts = 40

        if filetype == 'powerlaw':

            ns,r1,r2 = ast.separatrix_nsr(xybounds[0],xybounds[1],npoints=npts)
            r,ns1,ns2 = ast.separatrix_rns(xybounds[2],xybounds[3],npoints=npts)

            ax.plot(ns1,r,'black',linestyle='dashed')
            ax.plot(ns2,r,'black',linestyle='dashed')

            if fillste:
                ax.fill_betweenx(r,xybounds[0],ns1, facecolor='lightsalmon', alpha=0.2)
                p1 = Polygon([[0,0],[0,1]],color='lightsalmon',alpha=0.2)

                ax.fill_betweenx(r,ns1,ns2,facecolor='gold',alpha=0.3)
                p2 = Polygon([[0,0],[0,1]],color='gold',alpha=0.3)

                ax.fill_betweenx(r,ns2,xybounds[1],facecolor='skyblue',alpha=0.2)
                p3 = Polygon([[0,0],[0,1]],color='skyblue',alpha=0.2)

        elif filetype == 'slowroll':
            ste2,ste1 = ast.separatrix_epstwo_twoepsone(1e-10,xybounds[1],npoints=2*npts)

            if fillste:
                ax.fill_betweenx(ste1,ste2,xybounds[1], facecolor='lightsalmon', alpha=0.2)
                p1 = Polygon([[0,0],[0,1]],color='lightsalmon',alpha=0.2)

                ax.fill_between(ste2,ste1,xybounds[3],facecolor='gold',alpha=0.3)
                p2 = Polygon([[0,0],[0,1]],color='gold',alpha=0.3)

                ax.fill_betweenx([xybounds[2],xybounds[3]],xybounds[0],0.0,facecolor='skyblue',alpha=0.2)
                p3 = Polygon([[0,0],[0,1]],color='skyblue',alpha=0.2)

        else:
            raise Exception("create_figure: filetype unknown!")
        
        if fillste:
            ax.legend([p1,p2,p3], [r'$\epsilon_2 > 2 \epsilon_1$'
                                   , r'$0<\epsilon_2< 2 \epsilon_1$'
                                   , r'$\epsilon_2< 0$'] , loc='lower right', shadow=True
                      , fancybox=True, frameon=True, labelspacing=0
                      , fontsize=10)


    if contour or fillcontour:
        add_contour_plot(ax,contourfile[0],'cyan',contourfile[1],'royalblue',alpha=0.7
                         ,depth=2,fill=fillcontour,log=True)


    ax.set_xlabel(xylabel[0])
    ax.set_ylabel(xylabel[1])

    sc=add_model_scatter(ax,xycscat,scmap,marker='o',alpha=scalpha)
    add_model_annotate(ax,sannotate,xyannotate,xyptoffset,cstyle=('arc3,rad='+str(arc)))
    

    cb = fig.colorbar(sc,pad=0.03,aspect=30)
    cb.set_label(aio.paramdict['logEGeV'])

    
    
    if wrehbar:
        add_model_name(ax,model.getlabel()+r' $\mathrm{&}$ $\overline{w}_\mathrm{reh} =$'+ r''
                       +str(wrehbar))
    else:
        add_model_name(ax,model.getlabel()+r' $\mathrm{&}$ $\overline{w}_\mathrm{reh} = 0$')

    add_model_fixparams(ax,allfixed)    
    
    plt.savefig(figfile,bbox_inches='tight')
    plt.close()

