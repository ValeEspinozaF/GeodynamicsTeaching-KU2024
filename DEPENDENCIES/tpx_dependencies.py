# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 19:02:48 2023

@author: Valentina Espinoza
"""

import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable


def load_pickle_obj(file_path):
    with open(file_path, 'rb') as datafile:
        dict_data = pickle.load(datafile)
        
    return dict_data
    


def load_TPX(tpx_path):
    dict_data = load_pickle_obj(tpx_path)
    mdat = dict_data["mdat"][::-1]
    mlon = dict_data["mlon"]
    mlat = dict_data["mlat"]
    return mdat, mlon, mlat


def add_colorbar(fig, ax, im, label):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.1)
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(label, rotation=270, fontsize=11, labelpad=15)