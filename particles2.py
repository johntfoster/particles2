#!/usr/bin/env python

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
from lxml import etree
import os.path


class particle_realization():
    """
       Class to create a realization of densly packed circular particles (2D)
    """

    def __init__(self, width, height, particle_diameter, target_density=0.6,
                 driver_type='sine', sine_amp=1.0, sine_freq=1.0,
                 number_of_peridigm_nodes_across_particle_diameter=5):

        if height < sine_amp:
            print "There are no particles above the driver, please adjust the \
                   particle height or the sine wave amplitude"

        self.particle_diameter = particle_diameter
        self.particle_radius = particle_diameter / 2.0

        self.number_of_peridigm_nodes_across_particle_diameter = \
            number_of_peridigm_nodes_across_particle_diameter

        self.width = width
        self.height = height
        self.target_density = target_density
        self.driver_type = driver_type
        self.time = 0.0
        self.sin_amp = sine_amp
        self.sin_freq = sine_freq

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
        self.create_peridigm_walls()

    def __compute_total_area(self):

        full_area = self.width * self.height

        if self.driver_type == 'sine':
            driver_area, _ = scipy.integrate.quad(lambda x:
                                                  self.sin_amp * np.sin(self.sin_freq * x + np.arcsin(1.0))
                                                  + self.sin_amp, 0.0, self.width)
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

        particles = grid_pairs[grid_pairs[:,1] > (self.sin_amp * np.sin(self.sin_freq * grid_pairs[:,0] + np.arcsin(1.0)) + self.sin_amp + self.particle_diameter)]
        self.driver = grid_pairs[grid_pairs[:,1] <= (self.sin_amp * np.sin(self.sin_freq * grid_pairs[:,0] + np.arcsin(1.0)) + self.sin_amp + self.particle_diameter)]

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
              self.number_of_peridigm_nodes_across_particle_diameter / 2.0 )

        x = np.arange(0.0, self.width, dr)
        y = self.sin_amp * np.sin(self.sin_freq * x +
                                  np.arcsin(1.0)) + self.sin_amp

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

    def create_peridigm_walls(self):

        dr = (self.particle_diameter /
              self.number_of_peridigm_nodes_across_particle_diameter / 2.0)

        single_wall_grid = np.mgrid[0.0:5*dr:dr,
                                    0.0:(self.height):dr]

        wall_start_left_y = (self.sin_amp * np.sin(self.sin_freq * -dr/2.0 +
                             np.arcsin(1.0)) + self.sin_amp)

        wall_start_right_y = (self.sin_amp * np.sin(self.sin_freq *
                              -self.width + dr/2.0 +
                              np.arcsin(1.0)) + self.sin_amp)

        x_walls = np.array([single_wall_grid[0].ravel() - 5*dr - dr/2.0,
                            single_wall_grid[0].ravel() +
                            self.width + dr/2.0]).flatten()

        y_walls = np.array([single_wall_grid[1].ravel() + wall_start_left_y,
                            single_wall_grid[1].ravel() +
                            wall_start_right_y]).flatten()

        self.y_walls = y_walls[y_walls < self.height]
        self.x_walls = x_walls[y_walls < self.height]

    def plot_peridigm_nodes(self):


        plt.plot(self.x_driver, self.y_driver, 'ro', self.x, self.y, 'bo', 
                 self.x_walls, self.y_walls, 'go')
        plt.show()

    def print_peridigm_files(self, basename='peridigm'):

        f = open(basename+"_nodes.txt", 'w')
        g = open(basename+"_particles_nodeset.txt", 'w')
        h = open(basename+"_driver_nodeset.txt", 'w')
        i = open(basename+"_walls_nodeset.txt", 'w')

        node_id = 1
        f.write('#x y z block_id node_volume\n')
        for xy_loc in zip(self.x_driver, self.y_driver):

            f.write(str(xy_loc[0]) + " " + str(xy_loc[1]) + " 0.0 1 " +
                    str(self.node_volume) + "\n")
            h.write(str(node_id) + "\n")
            node_id += 1

        for xy_loc in zip(self.x_walls, self.y_walls):

            f.write(str(xy_loc[0]) + " " + str(xy_loc[1]) + " 0.0 2 " +
                    str(self.node_volume) + "\n")
            i.write(str(node_id) + "\n")
            node_id += 1

        for idx, particle_center in enumerate(zip(self.x, self.y)):

            particle_i = (self.particle_nodes +
                          [particle_center[0], particle_center[1]])

            for node in particle_i:
                f.write(str(node[0]) + " " + str(node[1]) + " 0.0 " +
                        str(idx+3) + " " + str(self.node_volume) + "\n")
                g.write(str(node_id) + "\n")
                node_id += 1

        f.close()
        g.close()
        h.close()
        i.close()

    def add_particle_output_to_xml(self, infile, freq=1, outfile=None, history_filename=None):

        if history_filename == None:
            history_basename = os.path.splitext(os.path.basename(infile))[0]
        else:
            history_basename = os.path.splitext(history_filename)[0]


        input_deck = etree.parse(infile)
        root = input_deck.getroot()

        #Create sublist for compute class
        compute_class = etree.SubElement(root, "ParameterList", name="Compute Class Parameters")

        #Create sublist for history output
        output2 = etree.SubElement(root, "ParameterList", name="Output2")
        etree.SubElement(output2, "Parameter", name="Output File Type",  type="string", value="ExodusII")
        etree.SubElement(output2, "Parameter", name="Output Format",  type="string",  value="BINARY")
        etree.SubElement(output2, "Parameter", name="Output Filename",  type="string",  value=history_basename)
        etree.SubElement(output2, "Parameter", name="Output Frequency",  type="int",  value=str(freq))
        etree.SubElement(output2, "Parameter", name="Parallel Write",  type="bool",  value="true")
        output_variables = etree.SubElement(output2, "ParameterList", name="Output Variables")
        
        count = 1
        for node in zip(self.x, self.y):
            velocity_tracker = etree.SubElement(compute_class, "ParameterList", name="Velocity of PD Node Nearest " + str(count))
            etree.SubElement(velocity_tracker, "Parameter", name="Compute Class", type="string", value="Nearest_Point_Data")
            etree.SubElement(velocity_tracker, "Parameter", name="X", type="double", value=str(node[0]))
            etree.SubElement(velocity_tracker, "Parameter", name="Y", type="double", value=str(node[1]))
            etree.SubElement(velocity_tracker, "Parameter", name="Z", type="double", value="0.0")
            etree.SubElement(velocity_tracker, "Parameter", name="Variable", type="string", value="Velocity")
            output_name = "Velocity_Particle_" + str(count)
            etree.SubElement(velocity_tracker, "Parameter", name="Output Label", type="string", value=output_name)
            etree.SubElement(output_variables, "Parameter", name=output_name, type="bool", value="true")
            count += 1

        #Write output to file
        if outfile == None:
            outname = os.path.splitext(os.path.basename(infile))[0] + "_full.xml"
        else:
            outname = outfile
        
        with open(outname, 'w') as f:
            input_deck.write(f, pretty_print=True)

real = particle_realization(1.0, 1.5, 0.02, target_density=0.6, sine_amp=0.25,
                            sine_freq=3.1415,number_of_peridigm_nodes_across_particle_diameter=4)
real.print_peridigm_files()
real.add_particle_output_to_xml('./Perturbation_textfile.xml', freq=100)

#real.print_lammps_datafile()
#real.plot_lammps_particles()
