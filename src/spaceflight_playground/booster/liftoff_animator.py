#!/usr/bin/python3

##
# @file liftoff_animator.py
# @author Paul Daum
# @brief Takes a trajectory of a liftoff and creates a nice animation of a spacecraft.
##


import sys, os
sys.path.append(os.path.realpath('../../../../'))
sys.path.append(os.getcwd())

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import numpy as np


##
# @class liftoff_animator
# @brief The animator class
##
class liftoff_animator:

    ##
    # @brief Initialization procedure
    # @param trajectory A liftoff_trajectory instance
    # @param figNum Number of the figure (default = 0)
    def __init__(self, trajectory, figNum = 0):
    
        # Trajectory to plot
        self.T = trajectory.total_time
        self.N = trajectory.num_samples
        self.params = trajectory.rocket.params
        self.xs = trajectory.states_initial_guess
        self.us = trajectory.us
        self.ds = trajectory.ds

        # Some rocket geometry parameters
        body_width = 5
        body_height = 35.0
        nose_height = 5.0
        nozzle_height = 5
        nozzle_start_width = 1.5
        nozzle_end_width = 3

        # Some reference points at zero angle
        self.com = [0,0]
        self.nozzle_joint = [0,-self.params['L']]
        self.zero_gimbal_ref_bottom = [0,-body_height]
        self.real_gimbal_ref_bottom = [0,-body_height+2]

        # Rocket vertices at zero angle
        self.nose_tip = [0, body_height - self.params['L'] + nose_height]
        self.body_topleft = [-body_width/2.0, body_height-self.params['L']]
        self.body_topright = [body_width/2.0, body_height-self.params['L']]
        self.body_bottomleft = [-body_width/2.0, -self.params['L']]
        self.body_bottomright = [body_width/2.0, -self.params['L']]
        self.nozzle_topleft = [-nozzle_start_width/2.0, -self.params['L']+2]
        self.nozzle_topright = [nozzle_start_width/2.0, -self.params['L']+2]
        self.nozzle_bottomleft = [-nozzle_end_width/2.0, -self.params['L']-nozzle_height]
        self.nozzle_bottomright = [nozzle_end_width/2.0, -self.params['L']-nozzle_height]

        # Some plot parameters
        self.margin = 40 
        self.rocket_color = 'b'
        self.ref_line_color = 'b'
        self.act_line_color = 'b'

        # Configure figure
        self.fig = plt.figure(figNum)
        self.ax = self.fig.add_subplot(111)
        self.ax.grid()
        self.ax.set_aspect('equal')
        plt.axis([
            float(self.xs[0,0])-self.margin, # xmin
            float(self.xs[0,0])+self.margin, # xmax
            float(self.xs[1,0])-self.margin, # ymin
            float(self.xs[1,0])+self.margin  # ymax
        ])
        
        # == STATIC PLOT ELEMENETS ==
        # Telemetry data box in upper right corner
        self.telemetry_box = patches.Rectangle(
            xy = (0.7,0.8),
            width = 0.3,
            height = 0.2,
            fill = True,
            transform = self.ax.transAxes,
            color = 'darkblue'
        )

        # Telemetry data labels
        self.ax.text(s = 'Speed: ', x = 0.72, y = 0.93, color = 'red', transform = self.ax.transAxes)
        self.ax.text(s = 'Altitude: ', x = 0.72, y = 0.85, color = 'red', transform = self.ax.transAxes)
        
        # Dashed vertical reference line
        self.angle_ref = patches.Polygon(
            xy = [[0.5, 0.1], [0.5, 0.9]], 
            linestyle = 'dashdot', 
            transform = self.ax.transAxes, 
            closed = None, 
            fill = None
        )

        # Add those elements to the plot
        self.ax.add_patch(self.telemetry_box)
        self.ax.add_patch(self.angle_ref)

        # == DYNAMIC PLOT ELEMENTS ==
        # Polygon for rocket body and nose
        self.body = patches.Polygon(
            xy = [self.body_bottomleft, 
                  self.body_topleft, 
                  self.nose_tip,
                  self.body_topright, 
                  self.body_bottomright],
            fill = True,
        )
        
        # Polygon for rocket nozzle
        self.nozzle = patches.Polygon(
            xy = [self.nozzle_topleft, self.nozzle_topright, self.nozzle_bottomright, self.nozzle_bottomleft],
            facecolor = 'red'
        )
        
        # Dashed reference line for zero-gimbal
        self.gimbal_ref = patches.Polygon(
            xy = [self.com, self.zero_gimbal_ref_bottom],
            linestyle = 'dashdot',
            closed = None,
            fill = None,
            edgecolor = 'g'
        )
        
        # Full line for actual gimbal angle
        self.gimbal_dir = patches.Polygon(
            xy = [self.nozzle_joint, self.real_gimbal_ref_bottom],
            linestyle = 'solid',
            closed = None,
            fill = None,
            edgecolor = 'b'
        )

        # Arc for gimbal angle
        self.gimbal_ang = patches.Arc(
            xy = self.nozzle_joint,
            width = 15,
            height = 15,
            angle = 270,
            theta1 = 0,
            theta2 = 180/np.pi * float(self.us[1,0]),
            edgecolor = 'b',
            linewidth = '2',
        )
        
        # Full line for actual gimbal angle
        self.gimbal_dir = patches.Polygon(
            xy = [self.nozzle_joint, self.real_gimbal_ref_bottom],
            linestyle = 'solid',
            closed = None,
            fill = None,
            edgecolor = 'b'
        )

        # Arc for total angle
        self.rocket_ang = patches.Arc(
            xy = self.com,
            width = 15,
            height = 15,
            angle = 90,
            theta1 = 0,
            theta2 = 180/np.pi * float(self.xs[1,0]),
            edgecolor = 'b',
            linewidth = '2'
        )
        

        # TODO: Rocket thrust

        # Add all patches to the plot
        self.ax.add_patch(self.body)
        self.ax.add_patch(self.nozzle)
        self.ax.add_patch(self.gimbal_ref)
        self.ax.add_patch(self.gimbal_dir)
        self.ax.add_patch(self.gimbal_ang)
        self.ax.add_patch(self.rocket_ang)

    ## 
    # @brief Updates all the dynamic objects in the plot. Is used by
    #        matplotlib.animation
    # @param i Index of the data
    ##
    def animation_main(self,i):
        # Get the transformation of the origin relative to the axis
        ts = self.ax.transData
        coords_com = ts.transform(self.com)

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

        # Set rocket and gimbal angle
        self.gimbal_ang.theta2 = 180.0/np.pi * float(self.us[1,i])
        self.rocket_ang.theta2 = 180.0/np.pi * float(self.xs[4,i])

        # Apply the transformations (rotations) to the plot elements
        self.body.set_transform(ts + t_rot_com)
        self.nozzle.set_transform(ts + t_rot_com + t_rot_base) 
        self.gimbal_ref.set_transform(ts + t_rot_com)
        self.gimbal_dir.set_transform(ts + t_rot_com + t_rot_base)
        self.gimbal_ang.set_transform(ts + t_rot_com)

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