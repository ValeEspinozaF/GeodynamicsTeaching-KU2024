# -*- coding: utf-8 -*-
"""
Created on Thu Dec 28 19:40:20 2023

@author: Valentina Espinoza
"""

# Load modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


from tpx_dependencies import load_TPX, add_colorbar
from map_dependencies import set_map
from figure_dependencies import create_figure


def add_topography(fig, ax,
                   mountain="ANDES",
                   colormap="terrain",
                   map_range=[-5, 5]):

    # Load grids 
    tp_path = "DATA/TPX_TP_%s.obj" %mountain
    mdat, mlon, mlat = load_TPX(tp_path)
           
    
    # Create the colormap
    if colormap == "terrain":
        color1 = plt.cm.gist_earth(np.linspace(0.4, 0.9, 128))
        color2 = plt.cm.Blues_r(np.linspace(0.2, 0.9, 128))
        color3 = np.vstack((color2, color1))
        cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', color3)
    
    else:
        # Check more colormaps in https://matplotlib.org/stable/users/explain/colors/colormaps.html
        cmap = colormap
    
    
    # Plot mdat surface colormap
    extent = mlon[0,:][0] - 1, mlon[0,:][-1] + 1, mlat[:,0][-1] - 1, mlat[:,0][0] + 1
    im = ax.imshow(mdat/1e3, origin='lower', extent=extent, cmap=cmap, 
                   vmin=map_range[0], vmax=map_range[1])
    
    
    # Add colorbar
    add_colorbar(fig, ax, im, label="Elevation (km)")
        
    
    # Finish map
    set_map(fig, ax)

        


if __name__ == "__main__":

    fig, ax = create_figure()
    add_topography(fig, ax,
                   mountain="ANDES")