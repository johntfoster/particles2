#!/usr/bin/env python

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt


class particle_realization():
    """
       Class to create a realization of densly packed circular particles (2D)
    """

    def __init__(self, width, height, particle_diameter, target_density=0.6,
                 driver_type='sine',
                 number_of_peridigm_nodes_across_particle_diameter=10):

        self.particle_diameter = particle_diameter
        self.particle_radius = particle_diameter / 2.0

        self.number_of_peridigm_nodes_across_particle_diameter = \
            number_of_peridigm_nodes_across_particle_diameter

        self.width = width
        self.height = height
        self.target_density = target_density
        self.driver_type = driver_type
        self.time = 0.0

        if target_density > 0.91:
            print "Warning, the theoretical maximum packing density \
                   in 2D is 91%"

        #Create arrays for the x and y positions of the points
        self.grid = np.mgrid[self.particle_radius:width:particle_diameter,
                             self.particle_radius:(height -
                                                   self.particle_diameter):
                             particle_diameter]

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

        self.create_peridigm_driver()

    def __compute_total_area(self):

        full_area = self.width * self.height

        if self.driver_type == 'sine':
            driver_area, _ = scipy.integrate.quad(lambda x:
                                                  np.sin(x + np.arcsin(1.0))
                                                  + 1.0, 0.0, self.width)
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


    def __discretize_single_particle(self):

        nr = self.number_of_peridigm_nodes_across_particle_diameter
        dr = self.particle_diameter/ nr / 2.0
        ds = dr

        particle_nodes = []

        for ir in range(1, nr+1):

            r = (ir - 0.5) * dr
            perim = 2.0 * np.pi * r
            nt = int(np.round(perim / ds))
            dt = 2.0 * np.pi / nt

            for it in range(1, nt+1):

                t = (it - 1.0) * dt
                particle_nodes += [[r * np.cos(t), r * np.sin(t)]]

        self.particle_nodes = np.array(particle_nodes)
        self.node_volume = (np.pi * self.particle_radius *
                            self.particle_radius /
                            len(self.particle_nodes))

    def create_peridigm_driver(self):

        #Following:
        #[http://stackoverflow.com/questions/19117660/how-to-generate-
        #equispaced-interpolating-values]
        #to get evenly spaced points along the curve

        dr = (self.particle_diameter /
              self.number_of_peridigm_nodes_across_particle_diameter)

        x = np.arange(0.0, self.width, dr)
        y = np.sin(x + np.arcsin(1.0)) + 1.0 

        xd = np.diff(x)
        yd = np.diff(y)
        dist = np.sqrt(xd * xd + yd * yd)
        u = np.cumsum(dist)
        u = np.hstack([[0], u])

        t = np.arange(0., u.max(), dr)
        x_temp = np.interp(t, u, x)
        y_temp = np.interp(t, u, y)

        self.x_driver = np.array([x_temp for i in range(5)]).flatten()
        self.y_driver = np.array([y_temp - i * dr
                                 for i in range(5)]).flatten()

    def plot_peridigm_nodes(self):

        self.create_peridigm_driver()

        plt.plot(self.x_driver, self.y_driver, 'ro', self.x, self.y, 'bo')
        plt.show()

    def print_peridigm_files(self, basename='peridigm'):

        f = open(basename+"_nodes.txt", 'w')
        g = open(basename+"_particles_nodeset.txt", 'w')
        h = open(basename+"_driver_nodeset.txt", 'w')

        node_id = 1
        for xy_loc in zip(self.x_driver, self.y_driver):

            f.write('#x y z block_id node_volume\n')
            f.write(str(xy_loc[0]) + " " + str(xy_loc[1]) + " 0.0 1 " +
                    str(self.node_volume) + "\n")
            h.write(str(node_id) + "\n")
            node_id += 1

        for idx, particle_center in enumerate(zip(self.x, self.y)):

            particle_i = (self.particle_nodes +
                          [particle_center[0], particle_center[1]])

            for node in particle_i:
                f.write(str(node[0]) + " " + str(node[1]) + " 0.0 " +
                        str(idx+2) + " " + str(self.node_volume) + "\n")
                g.write(str(node_id) + "\n")
                node_id += 1

        f.close()
        g.close()
        h.close()


real = particle_realization(5, 5, 0.2, target_density=0.6,
                            driver_type='sine')
real.print_peridigm_files()
real.plot_peridigm_nodes()

#real.print_lammps_datafile()
#real.plot_lammps_particles()
