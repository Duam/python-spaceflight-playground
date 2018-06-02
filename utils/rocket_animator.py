#!/usr/bin/python3

##
# @file rocket_animator.py
# @author Paul Daum
# @brief Takes a trajectory of:
#        Spacecraft angles,
#        Spacecraft masses,
#        Thrust magnitudes,
#        Thrust angles
#        and creates a nice animation of a spacecraft.
##


import sys, os
sys.path.append(os.path.realpath('../'))
sys.path.append(os.getcwd())

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

## 
# @class rocket_animator
# @brief The animator class
##
class rocket_animator:

    ##
    # @brief Initialization procedure
    # @param params Plotting parameters dictionary
    # @param figNum Number of the figure (default = 0)
    def __init__(self, params, figNum = 0):

        # Horizon length and sample number
        self.T = params['T']
        self.N = params['N']

        
