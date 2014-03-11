#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
import scipy.integrate
import scipy.spatial
import matplotlib.pyplot as plt
import ctypes


_particles = np.ctypeslib.load_library('_particles','.')
_particles.transfer_momentum.restype = None
_particles.wall_contact.restype = None
_particles.transfer_momentum.argtypes = [ctypes.c_void_p,
                                         ctypes.c_void_p,
                                         ctypes.c_void_p,
                                         ctypes.c_void_p,
                                         ctypes.c_int,
                                         ctypes.c_void_p,
                                         ctypes.c_void_p,
                                         ctypes.c_double]

_particles.wall_contact.argtypes = [ctypes.c_void_p,
                                    ctypes.c_void_p,
                                    ctypes.c_void_p,
                                    ctypes.c_void_p,
                                    ctypes.c_int,
                                    ctypes.c_double,
                                    ctypes.c_double,
                                    ctypes.c_int,
                                    ctypes.c_double]
                                      

class particle_realization():
    """
       Class to create a realization of densly packed circular particles (2D)
    """

    def __init__(self, width, height, particle_diameter, target_density=0.6, driver_type='sine'):

        self.particle_diameter = particle_diameter
        self.particle_radius = particle_diameter / 2.0

        self.width = width
        self.height = height
        self.target_density = target_density
        self.driver_type = driver_type
        self.time = 0.0

        if target_density > 0.91:
            print "Warning, the theoretical maximum packing density in 2D is 91%"

        #Create arrays for the x and y positions of the points
        self.grid = np.mgrid[self.particle_radius:width:particle_diameter,
                        self.particle_radius:(height-self.particle_radius):particle_diameter]

        #Shift every other row over to create hexagonal packing
        for idx, item in enumerate(self.grid[1]):
            if idx % 2 != 0:
                item += self.particle_radius

        self.x = self.grid[0].ravel()
        self.y = self.grid[1].ravel()

        if driver_type == 'sine':

            self.__remove_particles_randomly()

        particle_density = self.__compute_particle_density()

        print "Target particle density is: " + str(target_density)
        print "Actual particle density is: " + str(particle_density)

        self.__discretize_single_particle()


    def __compute_total_area(self):

        full_area = self.width * self.height

        if self.driver_type == 'sine':
            driver_area, _ = scipy.integrate.quad(lambda x: np.sin(x + np.arcsin(1.0))
                    + 1.0 + self.particle_diameter, 0.0, self.width)
            return full_area - driver_area
        else:
            return full_area

    def __compute_particle_area(self):

        return (np.pi * self.particle_diameter * 
                self.particle_diameter / 4.0) * len(self.x)

    def __compute_particle_density(self):

        return self.__compute_particle_area() / self.__compute_total_area()

    def __compute_number_of_particles_to_remove(self):

        return (len(self.x) - 
                len(self.x) * self.target_density / 
                self.__compute_particle_density())


    def __remove_particles_randomly(self):
        
        grid_pairs = np.array([self.grid[0].ravel(), self.grid[1].ravel()]).T

        particles = grid_pairs[grid_pairs[:,1] > (np.sin(grid_pairs[:,0] + np.arcsin(1.0)) + 1.0 + self.particle_diameter)]
        self.driver = grid_pairs[grid_pairs[:,1] <= (np.sin(grid_pairs[:,0] + np.arcsin(1.0)) + 1.0 + self.particle_diameter)]

        self.x = particles[:,0]

        number_of_particles_to_remove = self.__compute_number_of_particles_to_remove()

        np.random.shuffle(particles)

        self.x = particles[:-number_of_particles_to_remove,0]
        self.y = particles[:-number_of_particles_to_remove,1]

    def plot_lammps_particles(self):

        x_driver = self.driver[:,0]
        y_driver = self.driver[:,1]

        plt.plot(self.x, self.y, 'bo', x_driver, y_driver, 'ro')
        plt.show()

    def print_lammps_datafile(self,filename='particles.txt'):

        with open(filename, 'w') as f:

            num_atoms = len(self.x) + len(self.driver)
            f.write("#Lammps data file\n\n")
            f.write(str(num_atoms) + " atoms\n\n")
            f.write("2 atom types\n\n")
            f.write("0.0 " + str(self.width) + " xlo xhi\n" )
            f.write("0.0 " + str(self.height) + " ylo yhi\n" )
            f.write("0.0 0.0 zlo zhi\n\n" )
            f.write("Atoms\n\n")

            for idx, xy_loc in enumerate(self.driver):
                f.write(str(idx+1) + " 1 " + str(self.particle_diameter) +
                        " 1.0 " + str(xy_loc[0]) + " " + str(xy_loc[1]) +
                        " 0.0\n")
            for idx, xy_loc in enumerate(zip(self.x,self.y)):
                f.write(str(idx+len(self.driver)+1) + " 2 " +
                        str(self.particle_diameter) + " 1.0 " +
                        str(xy_loc[0]) + " " + str(xy_loc[1]) + " 0.0\n")


    def __search_for_contact_neighbors(self):

        grid_pairs = np.array([self.x, self.y]).T
         
        self.tree = scipy.spatial.cKDTree(grid_pairs)
        #neighbors_temp = self.tree.query_ball_point(grid_pairs, 5.0 * self.particle_diameter)
        _, neighbors = self.tree.query(grid_pairs, k=100, p=2, distance_upper_bound=5.0*self.particle_diameter)

        neighbors = np.delete(np.where(neighbors ==  self.tree.n, -1, neighbors),0,1)
        #Find the maximum length of any family, we will use this to recreate 
        #the families array such that it minimizes masked entries.
        self.neighbor_length_list = np.array((neighbors != -1).sum(axis=1), dtype=np.int32)
        #Recast the families array to be of minimum size possible
        self.neighbors = np.array(ma.masked_equal(neighbors, -1).compressed(), dtype=np.int32)

    def __wall_contact(self):

        if self.driver_type == 'sine':
            sine_wave_bool = 0
        else:
            sine_wave_bool = 1

        _particles.wall_contact(self.x.ctypes.data_as(ctypes.c_void_p),
                                self.y.ctypes.data_as(ctypes.c_void_p),
                                self.velocity_x.ctypes.data_as(
                                    ctypes.c_void_p),
                                self.velocity_y.ctypes.data_as(
                                    ctypes.c_void_p),
                                len(self.x), self.width, self.height,
                                sine_wave_bool, self.particle_radius)

    def initialize_velocities(self, min_velocity, max_velocity, scale_factor):

        self.velocity_x = ((scale_factor * max_velocity -
                           scale_factor * min_velocity) *
                           np.random.random_sample(len(self.x),)
                           + scale_factor * min_velocity)
        self.velocity_y = ((scale_factor * max_velocity -
                           scale_factor * min_velocity) *
                           np.random.random_sample(len(self.y),)
                           + scale_factor * min_velocity)

    def advance(self, dt):

        self.x += dt * self.velocity_x
        self.y += dt * self.velocity_y

    def transfer_momentum(self):
        """
           Calls C function to transfer momentum between particle contacts
        """

        _particles.transfer_momentum(
            self.x.ctypes.data_as(ctypes.c_void_p),
            self.y.ctypes.data_as(ctypes.c_void_p),
            self.velocity_x.ctypes.data_as(ctypes.c_void_p),
            self.velocity_y.ctypes.data_as(ctypes.c_void_p),
            len(self.x),
            self.neighbor_length_list.ctypes.data_as(ctypes.c_void_p),
            self.neighbors.ctypes.data_as(ctypes.c_void_p),
            self.particle_diameter)

        self.__wall_contact()

    def __discretize_single_particle(self, nr=10):

        dr = self.particle_radius / nr 
        ds = dr

        self.particle_nodes = []

        for ir in range(1, nr+1):

            r = (ir - 0.5) * dr
            perim = 2.0 * np.pi * r
            nt = np.rnd(perim / ds)
            dt = 2.0 * np.pi / nt

            for it in range(1, nt+1):

                t = (it - 1.0) * dt
                self.particle_nodes += [r * np.cos(t), r * np.sin(t)]

    def relax_particles(self, total_time):

        if self.time == 0.0:

            self.__search_for_contact_neighbors()

            self.initialize_velocities(-self.particle_diameter,
                                       self.particle_diameter, 1.0)

        while self.time < total_time:

            print("time: " + str(self.time))

            max_velocity = np.max(np.sqrt(self.velocity_x * self.velocity_x +
                                  self.velocity_y * self.velocity_y))
            dt_max = self.particle_radius / max_velocity / 2.0
            #dt_max = 1.e-6
            print("dt_max: " + str(dt_max))

            self.advance(dt_max)
            self.transfer_momentum()
            self.time = self.time + dt_max
            print self.time

    def animate_particle_motion(self, total_time=1.0, num_plots=10):

        plt.ion()
        data, = plt.plot(self.x, self.y, 'ro')
        plt.show()

        for time in np.arange(0, total_time, total_time/num_plots):
            real.relax_particles(time)
            data.set_xdata(self.x)
            data.set_ydata(self.y)
            plt.draw()


real = particle_realization(20, 10, 0.2, target_density=0.6,
                            driver_type='sine')

#real.print_lammps_datafile()
#real.plot_lammps_particles()
real.animate_particle_motion(total_time=0.001, num_plots=10000)
