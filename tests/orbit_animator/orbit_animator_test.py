#!/usr/bin/python3

## 
# @file orbit_animator_test.py
# @author Paul Daum
##

import sys, os
sys.path.append(os.path.realpath('../../'))
sys.path.append(os.getcwd())

import numpy as np

from models.kepler_orbit.kepler_orbit import kepler_orbit as orbit
from utils.xml_writer import read_from_xml
from utils.orbit_animator import orbit_animator


# Read the trajectory from the xml file
params, target_orbit, xs_out, us_out = read_from_xml('trajectory.xml')

orb_tar = orbit()
orb_tar.setOrbitalElements(
    e = np.array([target_orbit['e_x'], target_orbit['e_y'], 0]),
    h = np.array([0,0, target_orbit['h']])
)

print(xs_out[0,:])
print(xs_out[1,:])

#for k in range(params['N']):
#    xs_out[0,k] = xs_out[0,k] + 1737.5e3 

# Prepare parameter struct for the animator
anim_params = {
    'T': params['T'],
    'N': params['N'],
    'body_radius': 1737.5 * 1e3, #TODO add this info to the xml file
    'target_orbit': orb_tar,
    'isCartesian': True,
    'xPositions': xs_out[0,:],
    'yPositions': xs_out[1,:],
    'xVelocities': xs_out[2,:],
    'yVelocities': xs_out[3,:],
    'xForces': us_out[0,:],
    'yForces': us_out[1,:]
}

# Create an animator
anim = orbit_animator(anim_params)

# Run the animator
anim.run(10)
