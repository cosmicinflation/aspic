#A simple module to dump 2D contour lines that can be used with
#fooiplots from getdist samples

import numpy as np
import getdist.mcsamples as gm
import matplotlib.pyplot as plt



def contour_curves_2d(rootname,samples,paramx, paramy):
    
    density = samples.get2DDensityGridData(paramx,paramy, num_plot_contours=2,
                                       get_density=False, meanlikes=False)

    levels = np.flip(density.contours)
    
    cs = plt.contour(density.x,density.y,density.P,levels)


    nlevels = len(cs.allsegs)    
    
    for i in range(nlevels):
        csvertex = cs.allsegs[i][:]
        xycs = csvertex[0][:]
        xcs = np.array(xycs[:,0])
        ycs = np.array(xycs[:,1])

        todump = np.column_stack([xcs,ycs])
        filename = rootname + '_'+paramx+'_'+paramy+'_level_'+ str(i)+'.dat'
        np.savetxt(filename, todump)
    
        
