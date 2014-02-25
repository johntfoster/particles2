#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
import scipy.integrate
import scipy.spatial
import matplotlib.pyplot as plt


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

        self.__remove_particles_randomly()

        particle_density = self.__compute_particle_density()

        print("Target particle density is: " + str(target_density))
        print("Actual particle density is: " + str(particle_density))

        self.__search_for_contact_neighbors()

        self.initialize_velocities(-self.particle_diameter, self.particle_diameter)
        

    def __compute_total_area(self):

        full_area = self.width * self.height

        if self.driver_type == 'sine':
            driver_area, _ = scipy.integrate.quad(lambda x: np.sin(x + np.arcsin(1.0))
                    - 1.0, 0.0, self.height)
            return full_area - driver_area
        else:
            return full_area

    def __compute_particle_area(self):

        return (np.pi * self.particle_diameter / 4.0) * len(self.x)

    def __compute_particle_density(self):

        return self.__compute_particle_area() / self.__compute_total_area()

    def __compute_number_of_particles_to_remove(self):

        return (len(self.x) - 
                len(self.x) * self.target_density / 
                self.__compute_particle_density())

    def __remove_particles_randomly(self):

        number_of_particles_to_remove = self.__compute_number_of_particles_to_remove()

        grid_pairs = np.array([self.grid[0].ravel(), self.grid[1].ravel()]).T

        np.random.shuffle(grid_pairs)

        self.x = grid_pairs[:-number_of_particles_to_remove,0]
        self.y = grid_pairs[:-number_of_particles_to_remove,1]


    def __search_for_contact_neighbors(self):

        grid_pairs = np.array([self.x, self.y]).T
         
        self.tree = scipy.spatial.cKDTree(grid_pairs)

        neighbors_temp = self.tree.query_ball_point(grid_pairs, 5.0 * self.particle_diameter)

        length_of_neighbors = np.array([ len(item) for item in neighbors_temp ])

        max_neighbors = np.max(length_of_neighbors)

        self.neighbors = np.ones((len(neighbors_temp), max_neighbors), dtype=np.int) * -1

        for row in range(len(neighbors_temp)):
            self.neighbors[row,:length_of_neighbors[row]] = neighbors_temp[row]

        self.neighbors = ma.masked_equal(self.neighbors, -1) 
        self.neighbors.harden_mask()

    def __particles_in_contact(self):

        distances_x = ma.masked_array(self.x[self.neighbors] - self.x[:,None], 
                mask=self.neighbors.mask)
        distances_y = ma.masked_array(self.y[self.neighbors] - self.y[:,None], 
                mask=self.neighbors.mask)

        distances = np.sqrt(distances_x * distances_x + distances_y * distances_y)

        return (distances <= self.particle_diameter)

    def __wall_contact(self):

        case1 = (self.width - self.x) < self.particle_radius
        case2 = self.x < self.particle_radius
        case3 = (self.height - self.y) < self.particle_radius
        case4 = (self.y - np.sin(self.x + np.arcsin(1)) - 1) < self.particle_radius

        return (case1 | case2 | case3 | case4)


    def initialize_velocities(self,min_velocity, max_velocity):

        self.velocity_x = ((max_velocity - min_velocity) * 
                np.random.random_sample(len(self.x),) + min_velocity)
        self.velocity_y = ((max_velocity - min_velocity) * 
                np.random.random_sample(len(self.y),) + min_velocity)

    def advance(self, dt):

        self.x += dt * self.velocity_x
        self.y += dt * self.velocity_y

    def transfer_momentum(self):

        self.velocity_x = ma.masked_equal(np.where(self.__particles_in_contact(), self.velocity_x[self.neighbors], -1), -1).sum(axis=1)
        self.velocity_y = ma.masked_equal(np.where(self.__particles_in_contact(), self.velocity_y[self.neighbors], -1), -1).sum(axis=1)

        self.velocity_x = np.where(self.__wall_contact(), -self.velocity_x, self.velocity_x)
        self.velocity_y = np.where(self.__wall_contact(), -self.velocity_y, self.velocity_x)


    def relax_particles(self, total_time):

        while self.time < total_time:

            print("time: " + str(self.time))

            max_velocity = np.max([np.max(np.abs(self.velocity_x)), np.max(np.abs(self.velocity_y))])
            #dt_max = self.particle_radius / max_velocity / 2.0 / 10.0
            dt_max = 0.000000001
            print("dt_max: " + str(dt_max))

            self.advance(dt_max)
            self.transfer_momentum()

            self.time = self.time + dt_max
            print self.time

            


real = particle_realization(20,10,0.4,target_density=0.6)
plt.plot(real.x, real.y, 'ro')
plt.show()

real.relax_particles(0.1)
plt.plot(real.x, real.y, 'ro')
plt.show()
