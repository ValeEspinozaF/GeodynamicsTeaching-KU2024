# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 15:27:03 2024

@author: Valentina Espinoza
"""

# Public dependencies
import numpy as np
import pandas as pd


from figure_dependencies import create_figure, save_figure


def add_coastlines(figax=[],
                   line_color = "black",
                   line_opacity = 1.0,
                   line_width = 1.0,
                   save_plot = False):

    
    # Load coastlines
    file_path = "DATA/COASTLINE.txt"
    contour = np.loadtxt(file_path)

    
    # Create figure
    if figax:
        fig, ax = figax
    else:
        fig, ax = create_figure()
    
    
    ax.plot(contour[:,0], contour[:,1], 
            color=line_color, lw=line_width, alpha=line_opacity)
    
    
    # Save figure as png
    if save_plot:
        save_figure(save_plot)
        

    return fig, ax


add_coastlines(figax=[],
                   line_color = "black",
                   line_opacity = 1.0,
                   line_width = 1.0,
                   save_plot = False)