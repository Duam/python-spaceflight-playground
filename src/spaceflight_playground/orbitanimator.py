import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from spaceflight_playground.orbit_polar_model as orbit


class OrbitAnimator:
    def __init__(self, params, figNum = 0):
        """Takes a trajectory of position and velocities and tries to make a nice animation out of it.

        :param params: Plotting parameters dictionary
        :param figNum: Number of the figure
        """

        self.total_time = params['T']
        self.num_samples = params['N']
        self.planet_radius = params['body_radius']

        # Target orbit
        self.orbit_target = params['target_orbit']

        # Are the trajectories in cartesian coordinates?
        self.coordinates_are_cartesian = params['isCartesian']

        # If they're cartesian, take them as they are
        if self.coordinates_are_cartesian:
            self.xPositions = params['xPositions']
            self.yPositions = params['yPositions']
            self.xVelocities = params['xVelocities']
            self.yVelocities = params['yVelocities']
            self.xForces = params['xForces']
            self.yForces = params['yForces']

        # If they're polar, convert them
        else:
            rhos = params['rhos']
            thetas = params['thetas']
            rhoDots = params['rhoDots']
            thetaDots = params['thetaDots']
            rhoForces = params['rhoForces']
            thetaForces = params['thetaForces']

            # Stack them up into xs_pol, us_pol TODO
            # Convert into cartesian coordinates TODO

        # Configure figure
        margin = 2
        self.fig = plt.figure(figNum)
        self.ax = self.fig.add_subplot(
            111
            #xlim = (np.amin(self.xPositions)-margin, np.amax(self.xPositions)+margin),
            #ylim = (np.amin(self.yPositions)-margin, np.amax(self.yPositions)+margin)
        )
        self.ax.grid()
        self.ax.set_aspect('equal', adjustable='box')

        ## == STATIC PLOT ELEMENTS ==
        # Target orbit
        orbit_target_samples = self.orbit_target.discretize()
        orbit_target_plot, = self.ax.plot(
            orbit_target_samples[:,0], 
            orbit_target_samples[:,1], 
            '-', 
            color='red',
            lw=2
        )
        
        # Circle for celestial body
        celestial_body_plot = plt.Circle(
            (0,0), 
            self.planet_radius,
            color='grey', 
            lw=2, 
            fill=True
        )
        self.ax.add_patch(celestial_body_plot)


        ## == DYNAMIC PLOT ELEMENTS ==
        # Line plot for the trajectory
        self.trajectory = [[self.xPositions[0]],[self.yPositions[0]]]
        self.trajectory_plot, = self.ax.plot([], [],'o-', lw=2)

        # Osculating orbit
        self.orbit_osculating = orbit.kepler_orbit()
        self.orbit_osculating.fromCartesianState(
            self.xPositions[0],
            self.yPositions[0],
            self.xVelocities[0],
            self.yVelocities[0]
        )
        orbit_osculating_samples = self.orbit_osculating.discretize()
        self.orbit_osculating_plot, = self.ax.plot(
            orbit_osculating_samples[:,0],
            orbit_osculating_samples[:,1],
            '-',
            color='green',
            lw=2
        )


    def animation_main(self, i):
        """ Updates all the dynamic objects in the plot. Is used by matplotlib.animation.

        :param i: Index of the current frame
        :return: The updated dynamic objects.
        """

        # Line plot for trajectory
        self.trajectory[0].append(self.xPositions[i])
        self.trajectory[1].append(self.yPositions[i])
        self.trajectory_plot.set_data(self.trajectory[0], self.trajectory[1])
        
        self.orbit_osculating.fromCartesianState(
            self.xPositions[i],
            self.yPositions[i],
            self.xVelocities[i],
            self.yVelocities[i]
        )

        orbit_osculating_samples = self.orbit_osculating.discretize()
        self.orbit_osculating_plot.set_data(
            orbit_osculating_samples[:,0],
            orbit_osculating_samples[:,1]
        )

        return self.trajectory_plot#, self.orbit_osculating_plot


    def run(self, fps: float, saveFile: bool = False, fileName: str = 'anim.mp4') -> None:
        """Creates the animation and shows it.

        :param fps: Frames per second.
        :param saveFile: Saves the file if True.
        :param fileName: The filename.
        :return: None
        """

        # Create the animation object
        anim = animation.FuncAnimation(
            fig = self.fig,
            func = self.animation_main,
            frames = np.arange(1, self.num_samples),
            interval = 100/fps
        )

        # Save to file if needed
        if(saveFile):
            anim.save(fileName, fps)

        # Show animation
        plt.show()
