# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 17:19:01 2023

@author: Valentina Espinoza
"""

from tpx_dependencies import load_TPX, add_colorbar
from map_dependencies import set_map
from figure_dependencies import create_figure


def add_gravity(fig, ax,
                mountain="ANDES",
                colormap="bwr",
                map_range=[-600, 600]):

    # Load grids 
    tp_path = "DATA/TPX_GA_%s.obj" %mountain
    mdat, mlon, mlat = load_TPX(tp_path)       
    
    
    # Plot mdat surface colormap
    extent = mlon[0,:][0] - 1, mlon[0,:][-1] + 1, mlat[:,0][-1] - 1, mlat[:,0][0] + 1
    im = ax.imshow(mdat, origin='lower', extent=extent, cmap=colormap, 
                   vmin=map_range[0], vmax=map_range[1])
    
    
    # Add colorbar
    add_colorbar(fig, ax, im, label="Free-air gravity anomaly (mGal)")

    
    # Finish map    
    set_map(fig, ax)



if __name__ == "__main__":

    fig, ax = create_figure()
# =============================================================================
#     add_gravity(fig, ax,
#                 mountain="ANDES")
#     
# =============================================================================
    add_gravity(fig, ax,
                mountain="TIBET",
                colormap="gist_rainbow",
                map_range=[-200, 200])