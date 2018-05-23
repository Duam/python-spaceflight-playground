#!/usr/bin/python3

##
# @file conversion.py
# @author Paul Daum
##

import numpy as np

##
# @brief Converts cartesian coordinates to polar coordinates
# @param x The x coordinate in the cartesian frame
# @param y The y coordinate in the cartesian frame
# @return The polar coordinates [theta, rho]
##
def cart2pol (x, y):
    theta = np.arctan2(y, x)
    rho = np.hypot(x, y)
    return theta, rho

##
# @brief Converts polar coordinates to cartesian coordinates
# @param theta The angle coordinate in the polar frame
# @param rho The distance coordinate in the polar frame
# @return The cartesian coordinates [x, y]
##
def pol2cart (theta, rho):
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return x, y

## 
# @brief Converts a state vector in polar coordinates to cartesian coordinates
# @param x_pol The state vector in polar coordinates:
#           - rho: Distance from origin (m)
#           - theta: Angle, measured from horizontal plane in right-hand direction (rad)
#           - rhoDot. Radial velocity (m/s)
#           - thetaDot: Angular velocity (rad/s)
# @return The state vector in cartesian coordinates
#           - xPos: The x position (m)
#           - yPos: The y position (m)
#           - xVel: The x velocity (m/s)
#           - yVel: The y velocity (m/s)
##
def state_pol2cart (x_pol):

    # Get states
    rho = x_pol[0]      # Distance (Radius)
    theta = x_pol[1]    # Angle
    rhoDot = x_pol[2]   # Radial velocity
    thetaDot = x_pol[3] # Angular velocity

    # Compute position in cartesian frame
    xPos, yPos = pol2cart(theta, rho)

    # Expand position into a 3d vector (required for cross product)
    pos = np.array([xPos, yPos, 0])
    
    # Rewrite angular velocity as rotation around z-axis
    thetaDot = np.array([0, 0, thetaDot])

    # Compute radial and rotational velocities in cartesian frame
    vel_rad = np.array([rhoDot, 0, 0])
    vel_rot = np.cross(thetaDot, pos)

    # Add up velocities
    vel = vel_rad + vel_rot
    
    # Extract velocity components
    xVel = vel[0]
    yVel = vel[1]

    # Stack them into a vector
    x_cart = np.array([xPos, yPos, xVel, yVel])

    return x_cart


## 
# @brief Converts a state and control trajectory from polar coordinates
#        to cartesian coordinates
# @param xs_pol The state trajectory in polar coordinates:
#           - rho: Distance from origin (m)
#           - theta: Angle, measured from horizontal plane in right-hand direction (rad)
#           - rhoDot. Radial velocity (m/s)
#           - thetaDot: Angular velocity (rad/s)
# @param us_pol The control trajectory in polar coordinates
#           - u_rho The radial control force/acceleration
#           - u_theta The angular control force/acceleration
# @return All state vectors (stacked) in cartesian coordinates
#           - xPos: The x position (m)
#           - yPos: The y position (m)
#           - xVel: The x velocity (m/s)
#           - yVel: The y velocity (m/s)
# @return All control vectors (stacked) in cartesian coordinates
#           - u_x The control force/acceleration in x direction
#           - u_y The control force/acceleration in y direction
##        
def traj_pol2cart (xs_pol, us_pol):

    # Get the trajectory length
    N_x = xs_pol[:,0].size
    N_u = us_pol[:,0].size

    # Check if trajectories are compatible
    if (N_x != N_u):
        print("Error in traj_pol2cart: Trajectories are not equally long. Check your input.")
        return -1

    # Convert all state vectors to cartesian coordinates
    xs_cart = np.zeros([N_x, 4])
    for k in range(N_x):
        xs_cart = state_pol2cart(xs_pol[k,:])

    # Convert controls by using the corotating reference frame
    us_cart = us_pol
    for k in range(N_u):
        theta = xs_pol[k,1]
        rotation = np.array(
            [[np.cos(theta), -np.sin(theta)],
             [np.sin(theta),  np.cos(theta)]])
        us_cart[k,:] = np.dot(rotation, us_cart[k,:])

    # return
    return xs_cart, us_cart

##
# Execute this script to test the functions
##    
if __name__ == '__main__':
    print("NEEDS TESTING")