#!/usr/bin/python3

## 
# @file liftoff_trajectory.py
# @author Paul Daum
##

import sys, os
sys.path.append(os.path.realpath('../'))
sys.path.append(os.getcwd())

import xml.etree.cElementTree as etree

import numpy as np
import models.liftoff_model.liftoff_model


##
# @class liftoff_trajectory
# @brief Representation of a state and control 
#        trajectory of the liftoff_model class
##
class liftoff_trajectory:

    ##
    # @brief Initialization procedure
    # @param T Simulation time
    # @param N Number of samples
    # @param rocket liftoff_model instance (needed for parameters)
    ##
    def __init__(self, T, N, rocket = liftoff_model()):

        # Horizon length and sample number
        self.T = T
        self.N = N

        # Containers for states, controls, disturbances
        self.xs = np.zeros((rocket.nx, self.N+1))
        self.us = np.zeros((rocket.nu, self.N))
        self.ds = np.zeros((rocket.nd, self.N))

        # Set the initial state
        self.xs[:,0] = rocket.x0


    ##
    # @brief Setter for states
    # @param xs State trajectory with shape (N+1,nx)
    ##
    def setXs (self, xs):

        if(xs.shape != self.xs.shape)
            print("State trajectory shapes do not match.")
            return
        
        self.xs = xs
        return


    ##
    # @brief Setter for controls
    # @param us Control trajectory with shape (N,n)
    ##
    def setUs (self, us):

        if(us.shape != self.us.shape)
            print("Control trajectory shapes do not match.")
            return
        
        self.us = us
        return
    

    ##
    # @brief Setter for disturbances
    # @param us Control trajectory with shape (N,n)
    ##
    def setDs (self, ds):

        if(ds.shape != self.ds.shape)
            print("Disturbance trajectory shapes do not match.")
            return
        
        self.ds = ds
        return
    
    ##
    # @brief Adds a new state, control and disturbance
    # @param x State
    # @param u Control
    # @param d Disturbance
    # @param k Timestep
    ##
    def add (self, x, u, d, k):

        if(x.shape != self.xs[:,0].shape)
            print("State shapes do not match")
            return

        if(u.shape != self.us[:,0].shape)
            print("Control shapes do not match")
            return

        if(d.shape != self.ds[:,0].shape)
            print("Disturbance shapes do not match")
            return

        self.xs[:,k] = x
        self.us[:,k] = u
        self.ds[:,k] = d

        
    ## 
    # @brief Writes the trajectory to XML
    # @param filename The name of the xml file
    ##
    def writeXML(self, filename):
        
        # Create root and main branches
        root = etree.Element("Liftoff trajectory")
        params = etree.Element(root, 'Parameters')
        xs = etree.SubElement(root, 'xs')
        us = etree.SubElement(root, 'us')
        ds = etree.SubElement(root, 'ds')

        # Fill in parameters
        for key, value in self.params.items():
            etree.SubElement(params, key).text = str(value)

        # Fill in states
        for k in range(N+1):
            etree.SubElement(xs, 'x_'+str(k))
            for i in range(self.rocket.nx):
                etree.SubElement('x_'+str(k), self.rocket.x_keys[i]).text = str(x[k])

        # Fill in controls
        for k in range(N):
            etree.SubElement(us, 'u_'+str(k))
            for i in range(self.rocket.nu):
                etree.SubElement('u_'+str(k), self.rocket.u_keys[i]).text = str(u[k])

        # Fill in disturbances
        for k in range(N):
            etree.SubElement(xs, 'd_'+str(k))
            for i in range(self.rocket.nd):
                etree.SubElement('d_'+str(k), self.rocket.d_keys[i]).text = str(d[k])

        # Create tree and write to file
        tree = etree.ElementTree(root)
        tree.write(filename)


    # TODO: read from XML
    