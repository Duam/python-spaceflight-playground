#!/usr/bin/python3

##
# @file rocket_animator.py
# @author Paul Daum
# @brief Takes a trajectory of a liftoff and creates a nice animation of a spacecraft.
##


import sys, os
sys.path.append(os.path.realpath('../'))
sys.path.append(os.getcwd())

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import numpy as np

from models.liftoff_model.liftoff_trajectory import liftoff_trajectory

## 
# @class rocket_animator
# @brief The animator class
##
class rocket_animator:

    ##
    # @brief Initialization procedure
    # @param trajectory A liftoff_trajectory instance
    # @param figNum Number of the figure (default = 0)
    def __init__(self, trajectory, figNum = 0):
    
        # Trajectory to plot
        self.T = trajectory.T
        self.N = trajectory.N
        self.params = trajectory.rocket.params
        self.xs = trajectory.xs
        self.us = trajectory.us
        self.ds = trajectory.ds

        # Some rocket geometry parameters
        self.body_width = 5
        self.body_height = 35.0
        self.nose_height = 5.0
        self.nozzle_height = 5
        self.nozzle_start_width = 1.5
        self.nozzle_end_width = 3

        # Configure figure
        self.fig = plt.figure(figNum)
        self.ax = self.fig.add_subplot(111)
        self.ax.grid()
        self.ax.set_aspect('equal')
        plt.axis(
            [-40, 40, -40, 40]
        )
        
        # == STATIC PLOT ELEMENETS ==
        # Box in upper left corner for showing velocity

        # Dashed vertical reference line

        # == DYNAMIC PLOT ELEMENTS ==
        # Box for rocket body
        self.body = patches.Rectangle(
            xy = (-self.body_width/2.0, -self.params['L']),
            width = self.body_width,
            height = self.body_height,
            fill = True,
        )
        self.ax.add_patch(self.body)
    
        # Polygon for rocket nose
        nose_points = [
            [-self.body_width/2.0, self.body_height-self.params['L']], # Upper left body corner
            [0, self.body_height - self.params['L'] + self.nose_height], # Nose tip
            [self.body_width/2.0, self.body_height-self.params['L']] # Upper right body corner
        ]
        self.nose = patches.Polygon(
            xy = nose_points
        )
        self.ax.add_patch(self.nose)

        # Polygon for rocket nozzle
        nozzle_points = [
            [-self.nozzle_start_width/2.0, -self.params['L']+2], # Top left nozzle corner
            [-self.nozzle_end_width/2.0, -self.params['L']-self.nozzle_height], # Bottom left nozzle corner
            [self.nozzle_end_width/2.0, -self.params['L']-self.nozzle_height], # Bottom right nozzle corner
            [self.nozzle_start_width/2.0, -self.params['L']+2] # Top right nozzle corner
        ]
        self.nozzle = patches.Polygon(
            xy = nozzle_points,
            facecolor = 'red'
        )
        self.ax.add_patch(self.nozzle)

        # Dashed reference line for zero-gimbal
        self.gimbal_ref = patches.Polygon(
            xy = [[0,0],[0,-self.body_height]],
            linestyle = 'dashdot',
            closed = None,
            fill = None,
            edgecolor = 'g'
        )
        self.ax.add_patch(self.gimbal_ref)
        
        # Full line for actual gimbal angle
        # Arc for gimbal angle

    ## 
    # @brief Updates all the dynamic objects in the plot. Is used by
    #        matplotlib.animation
    # @param i Index of the data
    ##
    def animation_main(self,i):
        # Get the transformation of the origin relative to the axis
        ts = self.ax.transData
        coords_com = ts.transform([0,0])

        # Compute the transformation corresponding to a rotation
        # around the origin
        t_rot_com = mpl.transforms.Affine2D().rotate_around(
            x = coords_com[0], 
            y = coords_com[1], 
            theta = self.xs[4,i]
        )

        # Compute transformation for rotating around base
        coords_base = ts.transform([
            self.params['L'] * np.sin(self.xs[4,i]),
            -self.params['L'] * np.cos(-self.xs[4,i])
        ])
        t_rot_base = mpl.transforms.Affine2D().rotate_around(
            x = coords_base[0],
            y = coords_base[1],
            theta = self.us[1,i]
        )

        # Apply the transformations (rotations) to the plot elements
        self.nose.set_transform(ts + t_rot_com)
        self.body.set_transform(ts + t_rot_com)
        self.nozzle.set_transform(ts + t_rot_com + t_rot_base) 
        self.gimbal_ref.set_transform(ts + t_rot_com)

        return self.body, self.gimbal_ref, self.nose, self.nozzle

    ##
    # @brief Creates the animation and shows it
    ##
    def run(self, fps, saveFile = False, fileName = 'anim.mp4'):

        # Create the animation object
        anim = animation.FuncAnimation(
            fig = self.fig,
            func = self.animation_main,
            frames = np.arange(1,self.N),
            interval = 100/fps
        )

        # Save to file if needed
        if(saveFile):
            anim.save(fileName, fps)

        # Show animation
        plt.show()