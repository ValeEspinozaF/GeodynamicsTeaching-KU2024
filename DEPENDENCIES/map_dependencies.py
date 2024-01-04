# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 16:44:35 2023

@author: Valentina Espinoza
"""

import os
import numpy as np

from tpx_dependencies import load_TPX, add_colorbar


def set_map(fig, ax):
    
    setCartographic_AxisLabels(ax)
    ax.set_aspect('equal', 'box')
    ax.grid(alpha=0.3)



def add_coastlines(fig, ax,
                   line_color = "grey",
                   line_opacity = 1.0,
                   line_width = 1.0):
    
    # Check if ax has already data drawn in it
    has_drawn = False
    if any([ax.images, ax.lines, ax.collections]):
        has_drawn = True
        
    # Store original plot extent
    xmin, xmax= ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    
    # Load coastlines
    root_path = os.path.dirname(os.path.dirname(__file__))
    file_path = os.path.join(root_path, "DATA/COASTLINE.txt")
    contour = np.loadtxt(file_path)
    
    ax.plot(contour[:,0], contour[:,1], 
            color=line_color, lw=line_width, alpha=line_opacity)
    
    # Restore plot extent
    if has_drawn: 
        ax.set(xlim=(xmin, xmax), ylim=(ymin,ymax))



def add_plates(fig, ax,
               line_color = "black",
               line_opacity = 1.0,
               line_width = 1.2):
    
    # Check if ax has already data drawn in it
    has_drawn = False
    if any([ax.images, ax.lines, ax.collections]):
        has_drawn = True
    
    # Store original plot extent
    xmin, xmax= ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    
    # Load coastlines
    root_path = os.path.dirname(os.path.dirname(__file__))
    file_path = os.path.join(root_path, "DATA/PLATELINE.txt")
    contour = np.loadtxt(file_path)
    
    ax.plot(contour[:,0], contour[:,1], 
            color=line_color, lw=line_width, alpha=line_opacity)
    
    # Restore plot extent
    if has_drawn: 
        ax.set(xlim=(xmin, xmax), ylim=(ymin,ymax))



def add_topography(fig, ax,
                   mountain="ANDES",
                   colormap="jet",
                   map_range=[-5, 5]):

    # Load grids 
    root_path = os.path.dirname(os.path.dirname(__file__))
    file_path = os.path.join(root_path, "DATA/TPX_TP_%s.obj" %mountain)
    mdat, mlon, mlat = load_TPX(file_path)
    
    
    # Plot mdat surface colormap
    extent = mlon[0,:][0] - 1, mlon[0,:][-1] + 1, mlat[:,0][-1] - 1, mlat[:,0][0] + 1
    im = ax.imshow(mdat/1e3, origin='lower', extent=extent, cmap=colormap, 
                   vmin=map_range[0], vmax=map_range[1])
    
    
    # Add colorbar
    add_colorbar(fig, ax, im, label="Elevation (km)")
        
    
    # Finish map
    set_map(fig, ax)
    
    
    
def add_gravity(fig, ax,
                mountain="ANDES",
                colormap="bwr",
                map_range=[-600, 600]):

    # Load grids 
    root_path = os.path.dirname(os.path.dirname(__file__))
    file_path = os.path.join(root_path, "DATA/TPX_GA_%s.obj" %mountain)
    mdat, mlon, mlat = load_TPX(file_path)       
    
    
    # Plot mdat surface colormap
    extent = mlon[0,:][0] - 1, mlon[0,:][-1] + 1, mlat[:,0][-1] - 1, mlat[:,0][0] + 1
    im = ax.imshow(mdat, origin='lower', extent=extent, cmap=colormap, 
                   vmin=map_range[0], vmax=map_range[1])
    
    
    # Add colorbar
    add_colorbar(fig, ax, im, label="Free-air gravity anomaly (mGal)")

    
    # Finish map    
    set_map(fig, ax)


def setCartographic_AxisLabels(ax):
    
    # Ticks range
    ymin0, ymax0 = ax.get_ylim()
    xmin0, xmax0 = ax.get_xlim()
        
    
    # Ensure min starts at factor of 5 (if extent is not too small)
    if ymax0 - ymin0 > 5: 
        ymin, ymax = np.round( (ymin0-5)/5 ) * 5, ymax0
    else:
        ymin, ymax = np.floor(ymin0), np.ceil(ymax0)
        
        
    if xmax0 - xmin0 > 5: 
        xmin, xmax = np.round( (xmin0-5)/5 ) * 5, xmax0
    else:
        xmin, xmax = np.floor(xmin0), np.ceil(xmax0)
        
    
    # Set ticks
    yTickList = []
    xTickList = []
    step_sizes = [10, 5, 2, 1, 0.5]

    for step in step_sizes:
        yTickList = np.arange(ymin, ymax, step)
        if len(yTickList) >= 5:
            break
        
    for step in step_sizes:
        xTickList = np.arange(xmin, xmax, step)
        if len(xTickList) >= 5:
            break
    
    
    # Set labels
    xTickLabels = [""] * len(xTickList)
    for i in range(len(xTickList)):
        xtick = xTickList[i]
        
        if step >= 1:
            xTickLabels[i] = "%d$^\circ$" %np.abs(xtick)
        else:
            xTickLabels[i] = "%.1f$^\circ$" %np.abs(xtick)
        
        if xtick < 0:
            xTickLabels[i] += "W"
        elif xtick > 0:
            xTickLabels[i] += "E"
            
            
    yTickLabels = [""] * len(yTickList)
    for i in range(len(yTickList)):
        ytick = yTickList[i]
        
        if step >= 1:
            yTickLabels[i] = "%d$^\circ$" %np.abs(ytick)
        else:
            yTickLabels[i] = "%.1f$^\circ$" %np.abs(ytick)
        
        if ytick < 0:
            yTickLabels[i] += "S"
        elif ytick > 0:
            yTickLabels[i] += "N" 
    
    
    ax.set(
        yticks = yTickList,
        xticks = xTickList,
        yticklabels = yTickLabels,
        xticklabels = xTickLabels,
        ylim = (ymin0, ymax0),
        xlim = (xmin0, xmax0),
    )