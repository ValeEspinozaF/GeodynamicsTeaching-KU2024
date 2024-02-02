# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 15:47:28 2024

@author: Valentina Espinoza
"""

import numpy as np
import matplotlib.pyplot as plt



def elastic_bending_degree(elastic_thickness, density_crust=2700, density_mantle=3300):
    """
    This function quantifies the deflection degree of an elastic beam placed 
    on top of a viscous fluid. The maximum load is equivalent to the weight of 
    a 1-km-thick crustal block. Based on Turcotte & Schubert (2014), Equation 5.125.

    Parameters
    ----------
    elastic_thickness : float
        Plate's elastic thickness, expressed in meters.
    density_crust : float
        Density of crust, expressed in kg/m**3. Default is 2700 kg/m**3.
    density_mantle : float
        Density of mantle, expressed in kg/m**3. Default is 3300 kg/m**3.

    Returns
    -------
    degree : array
        2-column array where [0] contains the wavelength value (m), and 
        [1] contains the associated degree of deflection. Degree of deflection
        are values between 0 and 1, with 0 corresponding to no deflection, 
        while 1 corresponds to the maximum deflection. In other words, 0
        corresponds to the case where the load is supported fully through
        elasticity, while 1 corresponds to the case of full support through
        isostasy.

    """
    
    g = 9.8     # Gravity [m/s2]
    h0 = 1e3    # Load height [m]
    E = 0.5e11  # Young's modulus [Pa]
    nu = 0.2    # Poisson's ratio
    
    rho_c = density_crust
    rho_m = density_mantle
    Te = elastic_thickness

    wl = np.linspace(0, 3000e3, 1000)   # Wavelength [m]
    degree = np.zeros((len(wl), 2))
    degree[:, 0] = wl

    d1 = h0 / (rho_m / rho_c - 1 + np.zeros(len(wl)))
    
    D = (E * Te**3) / (12 * (1 - nu**2))        # Flexural rigidity
    A = D / (rho_c * g) * (2 * np.pi / wl)**4 
    d2 = h0 / (rho_m / rho_c - 1 + A)        # Deflection 

    degree[:, 1] = d2 / d1  # Degree of deflection

    return degree


def plot_elastic_bending_degree(elastic_thickness_values, density_crust=2700, density_mantle=3300):
    """
    Plots wavelength versus deflection degree of an elastic beam placed 
    on top of a viscous fluid. The maximum load is equivalent to the weight of 
    a 1-km-thick crustal block. Based on Turcotte & Schubert (2014), Equation 5.125.

    Parameters
    ----------
    elastic_thickness : float
        Plate's elastic thickness, expressed in meters.
    density_crust : float
        Density of crust, expressed in kg/m**3. Default is 2700 kg/m**3.
    density_mantle : float
        Density of mantle, expressed in kg/m**3. Default is 3300 kg/m**3.

    Returns
    -------
    fig, ax : matplotlib figure and axis objects
        Figure and axis objects containing the plot.
    """
    
    fig, ax = plt.subplots()
    
    
    ax.plot([0,3000], [1,1], '--', color="grey")
    ax.plot([0,3000], [0,0], '--', color="grey")
    ax.text(3000, 0, "no deflection", ha="right", va="bottom")
    ax.text(3000, 1, "max deflection", ha="right", va="top")
    
    title = ["Response of an elastic beam to a periodic load",
             "(max. push/pull equivalent to 1-km-thick crustal block)",
             rf"($\rho_C$ = {density_crust} $kg/m^3$, $\rho_M$ = {density_mantle} $kg/m^3$)",
             ]
    
    ax.set(title = "\n".join(title),
           xlabel = "Wavelength of periodic load [km]",
           ylabel = "Degree of vertical deflection",
           xlim = (0, 3000),
           )
    
    plt.grid(alpha=0.5)
    
    
    if isinstance(elastic_thickness_values, float):
        elastic_thickness_values = [elastic_thickness_values]
        
        
    for Te in elastic_thickness_values:

        g = 9.8     # Gravity [m/s2]
        h0 = 1e3    # Load height [m]
        E = 0.5e11  # Young's modulus [Pa]
        nu = 0.2    # Poisson's ratio
        
        rho_c = density_crust
        rho_m = density_mantle
        

        wl = np.linspace(0, 3000e3, 100000)
        elastic_bend_degree = np.zeros((len(wl), 2))
        elastic_bend_degree[:, 0] = wl

        d1 = h0 / (rho_m / rho_c - 1 + np.zeros(len(wl)))
        
        D = (E * Te**3) / (12 * (1 - nu**2))    # Flexural rigidity
        A = D / (rho_c * g) * (2 * np.pi / wl)**4   
        d2 = h0 / (rho_m / rho_c - 1 + A)

        elastic_bend_degree[:, 1] = d2 / d1
        
        
        
        # elastic thickness label
        x = wl / 1e3
        y = d2 / d1
        axis = [0, max(wl/1e3), 0, 1+0.5]
        dydx = np.diff(y) / np.diff(x)
        c = np.where(np.abs(dydx) == np.max(np.abs(dydx))) #
        indx = int(c[0])
        angle = float(np.arctan(dydx[indx] * np.diff(axis[0:2]) / np.diff(axis[2:4])) * 180 / np.pi)
        lbl = f"{Te/1e3} km"
        
        
        # plot data
        ax.plot(elastic_bend_degree[:,0]/1e3, elastic_bend_degree[:,1], '-')
        ax.text(np.mean(x[indx-1:indx]), np.mean(y[indx-1:indx]), lbl, 
                rotation=angle, ha="center", va="center", backgroundcolor="w")
        
        
    return fig, ax
        



if __name__ == "__main__":
    
    # Example usage:
    elastic_thickness = 35e3   # Average elastic thickness of the crust [m] 

    elastic_bend_degree = elastic_bending_degree(elastic_thickness)
    print("elastic_bend_degree: ", elastic_bend_degree)
    
    elastic_thickness_values = [35e3, 100e3]
    plot_elastic_bending_degree(elastic_thickness_values)
    
    
    
    
    
    