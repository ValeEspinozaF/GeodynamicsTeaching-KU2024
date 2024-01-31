# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 21:53:44 2024

@author: Valentina Espinoza
"""

import numpy as np
import pandas as pd
from math import sin, cos, asin, acos, atan2, sqrt, radians



def sph2cart(lon_deg, lat_deg, mag=1):
    """ Transforms spherical to cartesian coordinates
    
    Parameters
    ----------
    lon_deg : float, list or array
        Longitude in degrees.
    lat_deg : float, list or array
        Latitude in degrees.
    mag : float, list or array
        Angular magnitude. The default is 1.
    
    Returns
    ----------
    x_cart, y_cart, z_cart : floats or arrays
        Cartesian coordinates in same unit of measurement as mag.
    """
    
    lon_deg = np.asarray(lon_deg)
    lat_deg = np.asarray(lat_deg)

    if type(mag) in [int, float]:
        mag = np.full_like(lon_deg, mag)
    else:
        mag = np.asarray(mag)

    if type(lon_deg) != type(lat_deg):
        raise TypeError("Longitude and latitude input types are not the same.")

    if np.shape(lon_deg) != np.shape(lat_deg):
        raise ValueError("Longitude and latitude input sizes are not the same.")

    lon_rad = np.radians(lon_deg)
    lat_rad = np.radians(lat_deg)

    x_cart = mag * np.cos(lon_rad) * np.cos(lat_rad)
    y_cart = mag * np.cos(lat_rad) * np.sin(lon_rad)
    z_cart = mag * np.sin(lat_rad)

    return x_cart, y_cart, z_cart



def cart2sph(x_cart, y_cart, z_cart, out_fmt="deg"):
    """ Transforms cartesian to spherical coordinates.

    Parameters
    ----------
    x_cart, y_cart, z_cart : float, list or array
        DESCRIPTION.
    out_fmt : TYPE, optional
        Sets the format of the output coordinates (longitude and latitude). 
        The available options are "deg" for degrees and "rad" for radians. 
        The default is "deg".

    Returns
    -------
    lon, lat, mag : float or array
        Spherical coordinates with longitude and latitude in degrees or
        radians, and magnitude in the same unit of measurement as the input
        cartesian coordinates.
    """
    x_cart = np.asarray(x_cart)
    y_cart = np.asarray(y_cart)
    z_cart = np.asarray(z_cart)

    if type(x_cart) != type(y_cart) or type(y_cart) != type(z_cart):
        raise TypeError("X, Y, Z input types are not the same.")

    if np.shape(x_cart) != np.shape(y_cart) or np.shape(y_cart) != np.shape(z_cart):
        raise ValueError("X, Y, Z input sizes are not the same.")

    lon_rad = np.arctan2(y_cart, x_cart)
    lat_rad = np.arctan2(z_cart, np.sqrt(x_cart**2 + y_cart**2))
    mag = np.sqrt(x_cart**2 + y_cart**2 + z_cart**2)

    if out_fmt == "deg":
        lon_deg = np.degrees(lon_rad)
        lat_deg = np.degrees(lat_rad)
        return lon_deg, lat_deg, mag
    
    elif out_fmt == "rad":
        return lon_rad, lat_rad, mag
    
    else:
        raise ValueError("Invalid value for out_fmt. Use 'deg' or 'rad'.")



def ev_to_surfvel_eastnorth(ev_cart, pnt_lon, pnt_lat):
    """
    Calculates the east and north components of the surface velocity for
    a given spherical coordinate in the globe, and an Euler vector.

    Parameters
    ----------
    ev_cart : list or array
        Cartesian coordinates of the Euler vector (wx, wy, wz) in [deg/Myr].
    lat_deg : float or int
        Spherical latitude of a point on the Earth surface, in [degrees].
    lon_deg : float or int
        Spherical longitude of a point on the Earth surface, in [degrees].

    Returns
    -------
    v_east, v_north
        East and north components of the surface velocity in [cm/yr].

    """
    
    #Earth's radius
    Re = 6371000. 
    
    # Point in the earth in cartesian coordinates
    x, y, z = sph2cart(pnt_lon, pnt_lat, Re)  
    
    # Euler vector in radians/Myr
    ev_cart_array = np.atleast_2d(ev_cart)
    ev_cart_rad = ev_cart_array * (np.pi / 180) 
    
    # Surface velocity calculation
    sv_c_array = np.empty_like(ev_cart_array, dtype=np.float64)
    sv_c_array[:, 0] = (ev_cart_rad[:, 1] * z) - (ev_cart_rad[:, 2] * y)
    sv_c_array[:, 1] = (ev_cart_rad[:, 2] * x) - (ev_cart_rad[:, 0] * z)
    sv_c_array[:, 2] = (ev_cart_rad[:, 0] * y) - (ev_cart_rad[:, 1] * x)
    
    # [m/Myr] to [cm/yr]
    sv_c_array *= (1e2 / 1e6)
    
    lon_rad = np.radians(pnt_lon)
    lat_rad = np.radians(pnt_lat)

    vX = sv_c_array[:, 0]
    vY = sv_c_array[:, 1]
    vZ = sv_c_array[:, 2]

    R11 = -np.sin(lat_rad) * np.cos(lon_rad)
    R12 = -np.sin(lat_rad) * np.sin(lon_rad)
    R13 = np.cos(lat_rad)
    R21 = -np.sin(lon_rad)
    R22 = np.cos(lon_rad)
    R23 = 0.0

    v_north = R11 * vX + R12 * vY + R13 * vZ
    v_east = R21 * vX + R22 * vY + R23 * vZ

    return [float(v_east), float(v_north)] if len(ev_cart_array) == 1 else [v_east, v_north]



def eastnorth_to_azimuth(v_east, v_north):
    
    # Convert v_east to a NumPy array if it's a float
    if isinstance(v_east, (int, float)):
        v_east = np.array([v_east])

    # Convert v_north to a NumPy array if it's a float
    if isinstance(v_north, (int, float)):
        v_north = np.array([v_north])

    # Compute azimuth in degrees
    AT = np.degrees(np.arctan2(v_north, v_east))

    # Initialize azimuth array
    azimuth = np.empty_like(v_east)

    # Set values for v_east > 0
    r = np.where(v_east > 0)
    azimuth[r] = 90 - np.sign(v_north[r]) * np.abs(AT[r])

    # Set values for v_east <= 0
    r = np.where(v_east <= 0)
    azimuth[r] = -90 + np.sign(v_north[r]) * np.abs(AT[r])

    # Convert azimuth to float if it has only one element
    azimuth = float(azimuth) if len(azimuth) == 1 else azimuth

    return azimuth



def eastnorth_to_total(v_east, v_north):
    # Calculate total surface velocity
    
    v_total = np.sqrt(v_east**2 + v_north**2)
    return v_total


def ev_to_surfvel_total(ev_cart, pnt_lon, pnt_lat):
    """
    Calculates the surface velocity for a given spherical coordinate in 
    the globe, and an Euler vector.

    Parameters
    ----------
    ev_cart : list or array
        Cartesian coordinates of the Euler vector (wx, wy, wz) in [deg/Myr].
    lat_deg : float or int
        Spherical latitude of a point on the Earth surface, in [degrees].
    lon_deg : float or int
        Spherical longitude of a point on the Earth surface, in [degrees].

    Returns
    -------
    v_total
        Total surface velocity in [cm/yr].

    """
    
    v_east, v_north = ev_to_surfvel_eastnorth(ev_cart, pnt_lon, pnt_lat)
    v_total = eastnorth_to_total(v_east, v_north)
    return v_total


def ev_to_surfvel_azimuth(ev_cart, dir_azimuth, pnt_lon, pnt_lat):
    """
    Calculates the surface velocity along a particular direction
    for a given spherical coordinate in the globe, and an Euler vector.

    Parameters
    ----------
    ev_cart : list or array
        Cartesian coordinates of the Euler vector (wx, wy, wz) in [deg/Myr].
    dir_azimuth : float or int
        Direction azimuth, expressed in [degrees] clockwise from North.
    lat_deg : float or int
        Spherical latitude of a point on the Earth surface, in [degrees].
    lon_deg : float or int
        Spherical longitude of a point on the Earth surface, in [degrees].

    Returns
    -------
    dir_velocity
        Surface velocity along the given direction, in [cm/yr].

    """    


    v_east, v_north = ev_to_surfvel_eastnorth(ev_cart, pnt_lon, pnt_lat)
    vel_vector = np.array([v_east, v_north])
    
    dir_vector = np.array([ np.sin(np.radians(dir_azimuth)),  np.cos(np.radians(dir_azimuth)) ])
    dir_velocity = np.sum(vel_vector * dir_vector)
    
    return dir_velocity
    
    
    
def plot_surfvel_total(ax, torsvik_df, 
                       pnt_lon, pnt_lat, 
                       lineLabel = '',
                       lineParams = dict(),
                       ):

    time_list = torsvik_df["time1"].to_numpy()
    time_list = np.append(time_list, torsvik_df["time2"].iloc[-1])
    
    
    x_array = []
    y_array = []
    
    for i, stg_ev in torsvik_df.iterrows():    
        
        t1 = time_list[i]
        t2 = time_list[i + 1]
        x_array.extend([t1, t2])                
                
        ev_cart = sph2cart(stg_ev["lon"], stg_ev["lat"], stg_ev["angvel"])
        v_east, v_north = ev_to_surfvel_eastnorth(ev_cart, pnt_lon, pnt_lat)
        v_total = eastnorth_to_total(v_east, v_north)
        y_array.extend([v_total, v_total])

    
    
    # Plot line 
    lineDefault = {
        "color" : '0.2', 
        "linestyle" : '-', 
        "linewidth" : 1,
        "label" : lineLabel,
        "solid_capstyle" : "butt",
        "solid_joinstyle" : "miter",
    }

    # Update plot parameters   
    lineParams = {**lineDefault, **lineParams}
    ax.plot(x_array, y_array, **lineParams)    
    
    
    ax.set(xlabel = "Total velocity (cm/yr)",
           ylabel = "Time (Myr)")
    
    

def plot_surfvel_azimuth(ax, torsvik_df, 
                         pnt_lon, pnt_lat, 
                         lineLabel = '',
                         lineParams = dict(),
                         ):

    time_list = torsvik_df["time1"].to_numpy()
    time_list = np.append(time_list, torsvik_df["time2"].iloc[-1])
    
    
    x_array = []
    y_array = []
    
    for i, stg_ev in torsvik_df.iterrows():    
        
        t1 = time_list[i]
        t2 = time_list[i + 1]
        x_array.extend([t1, t2])                
                
        ev_cart = sph2cart(stg_ev["lon"], stg_ev["lat"], stg_ev["angvel"])
        v_east, v_north = ev_to_surfvel_eastnorth(ev_cart, pnt_lon, pnt_lat)
        v_total = eastnorth_to_azimuth(v_east, v_north)
        y_array.extend([v_total, v_total])

    
    
    # Plot line 
    lineDefault = {
        "color" : '0.2', 
        "linestyle" : '-', 
        "linewidth" : 1,
        "label" : lineLabel,
        "solid_capstyle" : "butt",
        "solid_joinstyle" : "miter",
    }

    # Update plot parameters   
    lineParams = {**lineDefault, **lineParams}
    ax.plot(x_array, y_array, **lineParams)    
    
    
    ax.set(xlabel = "Azimuth (degrees)",
           ylabel = "Time (Myr)")