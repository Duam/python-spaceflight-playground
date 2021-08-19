#!/usr/bin/python3

##
# @file rocket_animator.py
# @author Paul Daum
##

import sys, os
sys.path.append(os.path.realpath('../../'))
sys.path.append(os.getcwd())

from spaceflight_playground.booster import booster_model
from spaceflight_playground.booster import booster_trajectory
from spaceflight_playground.booster import booster_animation

import numpy as np


# Load a trajectory
spacecraft = booster_model()
trajectory = booster_trajectory(T=600.0, N=100)
trajectory.fromXML('orbit_animator_trajectory.xml')

trajectory.us[1,:] = np.pi/4 * np.ones((1,100))

# Create an animator
animator = booster_animation(trajectory=trajectory)

# Run the animator
animator.run(10)

'''
fig = plt.figure()
ax = fig.add_subplot(111)

body = patches.Rectangle(
    xy = (-1,-1),
    width = 2,
    height = 2,
    fill = True
)
ax.add_patch(body)

plt.show()'''