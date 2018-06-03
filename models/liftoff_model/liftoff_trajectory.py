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
import casadi as cas

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
    def __init__(self, T, N, rocket = None):

        if (rocket == None):
            rocket = liftoff_model()

        self.rocket = rocket

        # Horizon length and sample number
        self.T = T
        self.N = N

        # Containers for states, controls, disturbances
        self.xs = cas.DM.zeros((rocket.nx, self.N+1))
        self.us = cas.DM.zeros((rocket.nu, self.N))
        self.ds = cas.DM.zeros((rocket.nd, self.N))

        # Set the initial state
        self.xs[:,0] = rocket.x0


    ##
    # @brief Setter for states
    # @param xs State trajectory with shape (N+1,nx)
    ##
    def setXs (self, xs):

        if(xs.shape != self.xs.shape):
            print("State trajectory shapes do not match. xs.shape = " + str(xs.shape) + ", self.xs.shape = " + str(self.xs.shape))
            return
        
        self.xs = xs
        return


    ##
    # @brief Setter for controls
    # @param us Control trajectory with shape (N,n)
    ##
    def setUs (self, us):

        if(us.shape != self.us.shape):
            print("Control trajectory shapes do not match. us.shape = " + str(us.shape) + ", self.us.shape = " + str(self.us.shape))
            return
        
        self.us = us
        return
    

    ##
    # @brief Setter for disturbances
    # @param us Control trajectory with shape (N,n)
    ##
    def setDs (self, ds):

        if(ds.shape != self.ds.shape):
            print("Disturbance trajectory shapes do not match. ds.shape = " + str(ds.shape) + ", self.ds.shape = " + str(self.ds.shape))
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

        if(x.shape != self.xs[:,0].shape):
            print("State shapes do not match")
            return

        if(u.shape != self.us[:,0].shape):
            print("Control shapes do not match")
            return

        if(d.shape != self.ds[:,0].shape):
            print("Disturbance shapes do not match")
            return

        self.xs[:,k] = x
        self.us[:,k] = u
        self.ds[:,k] = d

        
    ## 
    # @brief Writes the trajectory to XML
    # @param filename The name of the xml file
    ##
    def toXML(self, filename):
        
        # Create root and main branches
        root = etree.Element("Liftoff trajectory")
        params = etree.SubElement(root, 'Parameters')
        xs = etree.SubElement(root, 'xs')
        us = etree.SubElement(root, 'us')
        ds = etree.SubElement(root, 'ds')

        # Fill in parameters
        for key, value in self.rocket.params.items():
            etree.SubElement(params, key).text = str(value)

        # Fill in states
        for k in range(self.N+1):
            x = {}
            for i in range(self.rocket.nx):
                x[self.rocket.x_keys[i]] = str(self.xs[i,k])

            etree.SubElement(xs, 'x', k = str(k), **x)
        
        # Fill in controls
        for k in range(self.N):
            u = {}
            for i in range(self.rocket.nu):
                u[self.rocket.u_keys[i]] = str(self.us[i,k])
            
            etree.SubElement(xs, 'u', k = str(k), **u)
        
        # Fill in disturbances
        for k in range(self.N):
            d = {}
            for i in range(self.rocket.nd):
                d[self.rocket.d_keys[i]] = str(self.ds[i,k])

            etree.SubElement(xs, 'd', k = str(k), **d)
        
        # Create tree and write to file
        tree = etree.ElementTree(root)
        tree.write(filename)


    # TODO: read from XML
    ##
    # @brief Reads the trajectory from an XML file
    # @param filename The name of the xml file
    ##
    def fromXML(self, filename):

        # Grab the root and main branches of the element tree
        root = etree.parse(filename).getroot()
        params = root.find('Parameters')
        xs = root.find('xs')
        us = root.find('us')
        ds = root.find('ds')

        # Grab parameters
        for key, value in params.items():
            self.params[key] = float(value)

            # Special case: N is an integer
            if (key == 'N'):
                self.params[key] = int(value)

        # Grab states
        for k in range(N+1):
            for i in range(self.rocket.nx):
                xs[i,k] = float(xs.find(self.rocket.x_keys[i]))

        # Grab controls
        for k in range(N):
            for i in range(self.rocket.nu):
                us[i,k] = float(us.find(self.rocket.u_keys[i]))

        # Grab disturbances
        for k in range(N):
            for i in range(self.rocket.nd):
                ds[i,k] = float(ds.find(self.rocket.d_keys[i]))
