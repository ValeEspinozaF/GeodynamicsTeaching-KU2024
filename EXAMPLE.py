# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 16:25:19 2024

@author: Valentina Espinoza
"""

# Load dependencies module
import os
import sys
dep_path = os.path.join(os. getcwd(), "DEPENDENCIES")
if not dep_path in sys.path:
    sys.path.append(dep_path)
    

from map_dependencies import add_plates, add_coastlines, add_gravity, add_topography
from figure_dependencies import create_figure

fig1, ax1 = create_figure()
add_gravity(fig1, ax1,
            mountain="TIBET",
            colormap="gist_rainbow",
            map_range=[-200, 200])

add_coastlines(fig1, ax1)


# =============================================================================
# 
# fig2, ax2 = create_figure()
# add_topography(fig2, ax2,
#                mountain="TIBET",
#                colormap="jet",
#                map_range=[0, 7])
# 
# 
# add_plates(fig2, ax2)
# =============================================================================
