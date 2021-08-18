import numpy as np
from spaceflight_playground.models.kepler_orbit import KeplerOrbit
from spaceflight_playground.utils import read_from_xml
from spaceflight_playground.orbiter.orbiter_animation import OrbitAnimator


# Read the trajectory from the xml file
params, target_orbit, xs_out, us_out = read_from_xml('orbit_animator_trajectory.xml')

orb_tar = KeplerOrbit()
orb_tar.setOrbitalElements(
    eccentricity=np.array([target_orbit['e_x'], target_orbit['e_y'], 0]),
    specific_angular_momentum=np.array([0, 0, target_orbit['h']])
)

print(xs_out[0, :])
print(xs_out[1, :])

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
anim = OrbitAnimator(anim_params)

# Run the animator
anim.run(10)


##
# Execute this script to test the animator
##
if __name__ == '__main__':


    """ From another file """
    # Create a trajectory
    T = 10.0
    N = 10
    R = 5

    xPoses = [0,1,2,3,4,5,6,7,8,9]
    yPoses = [9,8,7,6,5,4,3,2,1,0]
    xVelos = [1,1,1,1,1,1,1,1,1,1]
    yVelos = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]
    xForces = [0,0,0,0,0,0,0,0,0,0]
    yForces = [0,0,0,0,0,0,0,0,0,0]

    # Create a target orbit
    orbit_target = KeplerOrbit.kepler_orbit()
    orbit_target.fromEllipseParams(0.3, 0.0, 3)

    # Save trajectory in a dictionary
    params = {}
    params['T'] = T
    params['N'] = N
    params['body_radius'] = R
    params['target_orbit'] = orbit_target
    params['isCartesian'] = True
    params['xPositions'] = xPoses
    params['yPositions'] = yPoses
    params['xVelocities'] = xVelos
    params['yVelocities'] = yVelos
    params['xForces'] = xForces
    params['yForces'] = yForces

    # Create animator
    animator = OrbitAnimator(params)

    # Run animator
    animator.run(fps=1)

    # TODO load pre-optimized trajectory from xml