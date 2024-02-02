# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 16:57:07 2023

@author: Valentina Espinoza
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def create_figure(figsize=(5,5), dpi=100):
    """
    Create an empty figure and axes.

    Parameters
    ----------
    figsize : tuple, optional
        Size of the figure in (width, height). Default is (5,5).
    dpi : int, optional
        Dots per inches. Determines the resolution (how many pixels) 
        the figure comprises. Default is 100.

    Returns
    -------
    fig, ax : matplotlib figure and axis objects
        Objects necessary to modify and re-plot the figure.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    return fig, ax


def save_figure(fig_path, dpi=360):   
    """
    Save the figure to the specified path.
    """ 
    plt.savefig(fig_path, bbox_inches='tight', dpi=dpi)