# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 16:44:35 2023

@author: Valentina Espinoza
"""

import os
import numpy as np
from math import sin, cos, atan2, sqrt, radians

from tpx_dependencies import load_TPX, add_colorbar
from figure_dependencies import create_figure

def add_residual_bathymetry_hawaii(fig, ax,
                                   colormap="jet",
                                   map_range=[-2500, 2500]):
    """
    Add the residual bathymetry of the Hawaii region to a map.

    Parameters
    ----------
    fig : matplotlib figure object
        
    ax : matplotlib axis object
    
    colormap : str, optional
        Colormap keyword (accessible via matplotlib.colormaps). 
        Default is "jet".
    map_range : list, optional
        Minimum and maximum value for the colormap, expressed in 
        meters. Default is [-2500, 2500].

    Returns
    -------
    image : matplotlib image object
        Object necessary to extract information on the drawn topography.
    """

    # Load grid
    root_path = os.path.dirname(os.path.dirname(__file__))
    file_path = os.path.join(root_path, "DATA/RESIDUAL_BATHYMETRY_HAWAII.txt")
    grid_data = np.loadtxt(file_path)

    
    # Plot mdat surface colormap
    extent=[-170, -147, 12, 30]
    image = ax.imshow(grid_data, origin='upper', extent=extent, cmap=colormap, 
                      vmin=map_range[0], vmax=map_range[1])
    
    
    # Add colorbar
    add_colorbar(fig, ax, image, label="Residual bathymetry (m)")
        
    
    # Finish map
    #set_map(fig, ax)
    
    return image

def basic_map(figsize=(5,5), dpi=100,
              xlim=(-180, 180),
              ylim=(-90, 90)):
    """
    Create a basic map with coastlines and plate boundaries.

    Parameters
    ----------
    figsize : tuple, optional
        Size of the figure in (width, height). Default is (5,5).
    dpi : int, optional
        Dots per inches. Determines the resolution (how many pixels) 
        the figure comprises. Default is 100.
    xlim : tuple, optional
        Range of longitudes to be displayed in the figure. Default is (-180, 180).
    ylim : tuple, optional
        Range of latitudes to be displayed in the figure. Default is (-90, 90).

    Returns
    -------
    fig, ax : matplotlib figure and axis objects
        Objects necessary to modify and re-plot the map.
    """
    
    fig, ax = create_figure(figsize, dpi)
    add_plates(fig, ax)
    add_coastlines(fig, ax)
    ax.set(xlim=xlim, ylim=ylim)
    set_map(fig, ax)
    
    return fig, ax
    
    

def set_map(fig, ax):
    
    setCartographic_AxisLabels(ax)
    ax.set_aspect('equal', 'box')
    ax.grid(alpha=0.3)



def add_coastlines(fig, ax,
                   line_color = "grey",
                   line_opacity = 1.0,
                   line_width = 1.0):
    """
    Add coastlines to a map.

    Parameters
    ----------
    fig : matplotlib figure object
        
    ax : matplotlib axis object
        
    line_color : str, optional
        Color of the line. Default is "grey".
    line_opacity : float, optional
        Opacity of the line (values from 0 to 1, 1 is no transparency). 
        Default is 1.0.
    line_width : float, optional
        Width of the line. Default is 1.0.
    """
    
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
    """
    Add plate boundaries to a map.

    Parameters
    ----------
    fig : matplotlib figure object
        
    ax : matplotlib axis object
        
    line_color : str, optional
        Color of the line. Default is "grey".
    line_opacity : float, optional
        Opacity of the line (values from 0 to 1, 1 is no transparency). 
        Default is 1.0.
    line_width : float, optional
        Width of the line. Default is 1.0.
    """
    
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
                   mountain,
                   colormap="jet",
                   map_range=[-5, 5]):
    """
    Add topography to a map.

    Parameters
    ----------
    fig : matplotlib figure object
        
    ax : matplotlib axis object
    
    mountain : str, optional
        Keyword for the orogen. Available options are "ANDES", "CENTRAL_ROCKIES",
        "NORTH_ROCKIES", "SOUTH_ROCKIES", "TIBET" and "ZAGROS".
    colormap : str, optional
        Colormap keyword (accessible via matplotlib.colormaps). 
        Default is "jet".
    map_range : list, optional
        Minimum and maximum value for the colormap, expressed in 
        kilometers. Default is [-5, 5].

    Returns
    -------
    image : matplotlib image object
        Object necessary to extract information on the drawn topography.
    """

    # Load grids 
    root_path = os.path.dirname(os.path.dirname(__file__))
    file_path = os.path.join(root_path, "DATA/TPX_TP_%s.obj" %mountain)
    mdat, mlon, mlat = load_TPX(file_path)
    
    
    # Plot mdat surface colormap
    extent = mlon[0,:][0] - 1, mlon[0,:][-1] + 1, mlat[:,0][-1] - 1, mlat[:,0][0] + 1
    image = ax.imshow(mdat/1e3, origin='lower', extent=extent, cmap=colormap, 
                      vmin=map_range[0], vmax=map_range[1])
    
    
    # Add colorbar
    add_colorbar(fig, ax, image, label="Elevation (km)")
        
    
    # Finish map
    set_map(fig, ax)
    
    return image
    
    
    
def add_gravity(fig, ax,
                mountain="ANDES",
                colormap="bwr",
                map_range=[-600, 600]):
    """
    Add Free-air gravity anomaly to a map.

    Parameters
    ----------
    fig : matplotlib figure object
        
    ax : matplotlib axis object
    
    mountain : str, optional
        Keyword for the orogen. Available options are "ANDES", "CENTRAL_ROCKIES",
        "NORTH_ROCKIES", "SOUTH_ROCKIES", "TIBET" and "ZAGROS".
    colormap : str, optional
        Colormap keyword (accessible via matplotlib.colormaps). 
        Default is "jet".
    map_range : list, optional
        Minimum and maximum value for the colormap, expressed in mGal. 
        Default is [-5, 5].

    Returns
    -------
    image : matplotlib image object
        Object necessary to extract information on the drawn gravity anomaly.
    """

    # Load grids 
    root_path = os.path.dirname(os.path.dirname(__file__))
    file_path = os.path.join(root_path, "DATA/TPX_GA_%s.obj" %mountain)
    mdat, mlon, mlat = load_TPX(file_path)       
    
    
    # Plot mdat surface colormap
    extent = mlon[0,:][0] - 1, mlon[0,:][-1] + 1, mlat[:,0][-1] - 1, mlat[:,0][0] + 1
    image = ax.imshow(mdat, origin='lower', extent=extent, cmap=colormap, 
                      vmin=map_range[0], vmax=map_range[1])
    
    
    # Add colorbar
    add_colorbar(fig, ax, image, label="Free-air gravity anomaly (mGal)")

    
    # Finish map    
    set_map(fig, ax)
    
    return image



def geodesic_distance(point1_lon, point1_lat, 
                      point2_lon, point2_lat,
                      radius = 6371e3):
    
    """ Calculates the geodesic distance between two point on a sphere
    of given radius. """
    
    
    # Turn input coordinates from sph to radians
    lat1, lon1 = radians(point1_lat), radians(point1_lon)
    lat2, lon2 = radians(point2_lat), radians(point2_lon)

    a = cos(lat2)*sin(abs(lon2 - lon1))
    b = cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(abs(lon2 - lon1))
    c = sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(abs(lon2 - lon1))
    
    # Geodetic distance in meters
    distance = radius * atan2( sqrt(a**2 + b**2), c )
    
    return distance



def shifted_point(point1_lon, point1_lat, azimuth, distance=15e-4):
    """ Given a point on a sphere, generate another point located
    along a line defined by an azimuth angle in degrees. """
    
    # Shifted point in radians
    point2_lon = np.radians(point1_lon) + (np.radians(distance) * np.sin(np.radians(azimuth)))
    point2_lat = np.radians(point1_lat) + (np.radians(distance) * np.cos(np.radians(azimuth))) 
    
    # Shifted point in degrees
    point2_lat = np.degrees(point2_lat)
    point2_lon = np.degrees(point2_lon)
    
    return point2_lon, point2_lat


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
