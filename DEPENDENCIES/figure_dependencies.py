# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 16:57:07 2023

@author: Valentina Espinoza
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def create_figure(figsize=(10,7), dpi=360):
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    return fig, ax


def save_figure(fig_path, dpi=360):    
    plt.savefig(fig_path, bbox_inches='tight', dpi=dpi)