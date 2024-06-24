"""
Author: Blythe Fernandes
UUN number: s1967975

 CompMod Ex2: Particle3D, a class to describe point particles in 3D space

 An instance describes a particle in Euclidean 3D space: 
 velocity and position are [3] arrays

 Includes time integrator methods('update_pos' , 'update_pos_2nd' , 
 'update_vel'), string method (using '__str__'), methods for kinetic energy 
 and momentum.
 
 The static methods in this class are mentioned and described below under the 
 class Particle3D.
"""

import numpy as np

class Particle3D(object):
    """
    Class to describe point-particles in 3D space

        Properties:
    label : name of the particle
    mass  : mass of the particle
    pos   : position of the particle
    vel   : velocity of the particle

        Methods:
    __init__       - Initialises a particle properties in 3D space
    __str__        - Produces a string storing the particle properties
    kinetic_e      - computes the kinetic energy
    momentum       - computes the linear momentum
    update_pos     - updates the position to 1st order
    update_pos_2nd - updates the position to 2nd order
    update_vel     - updates the velocity

        Static Methods:
    new_p3d      - initializes a P3D instance from a file handle
    sys_kinetic  - computes total Kinetic Energy of a p3d list
    com_velocity - computes total mass and CoM (Center of Mass) velocity of a 
                   p3d list
    """
     
    def __init__(self, label, mass, pos, vel):
        """
        Initialises a particle in 3D space

        :param label: String,the name of the particle
        :param mass: float, mass of the particle
        :param pos: [3] float (numpy)array, position of the particle
        :param vel: [3] float (numpy)array, velocity of the particle
        """
        #Defining properties for this class
        self.label = str(label)
        self.mass = float(mass)
        
        #pos = float(input("Enter the x position of the particle: ")),float(input("Enter the y position of the particle: ")),float(input("Enter the z position of the particle: "))
        self.pos = np.array(pos, float) #Changing position vector into array float using numpy array
        
        #vel = float(input("Enter the x position of the particle: ")),float(input("Enter the y position of the particle: ")),float(input("Enter the z position of the particle: "))
        self.vel = np.array(vel, float) #Changing velocity vector into array float using numpy array
        

    def __str__(self):
        """
        XYZ-compliant string.
        
        :return (string) format: <label>    <x>  <y>  <z>
        """
        return (str(self.label) + "    " + str(self.pos[0]) + "   " + str(self.pos[1]) + "   " + str(self.pos[2]))


    def kinetic_e(self):
        """
        Returns the kinetic energy of a Particle3D instance

        :return ke: float, 1/2  * m * v**2
        """
        
        ke = 0.5 * self.mass * (np.square(self.vel[0]) + np.square(self.vel[1]) + np.square(self.vel[2]))
        #ke = np.dot((0.5 * self.mass), np.square(self.vel))
        return ke


    def momentum(self):
        """
        Returns the momentum of the Particles
        
        :return p: vector, m * vel(vector)
        """
        
        p = np.multiply(self.vel, self.mass) #Using a numpy function
        return p


    def update_pos(self, dt):
        """
        Updates the position vector(pos) to the 1st order that is within an 
        infinitesimal timestep for the change in the particles position.
        
        :param dt: float, timestep which is the change in time
        
        This function does not return anything
        """
        
        self.pos += np.multiply(self.vel, float(dt))


    def update_pos_2nd(self, dt, f):
        """
        Updates the position vector(pos) to the 2nd order that is within an 
        infinitesimal timestep and a force acting on the particle.
        
        :param dt: float, timestep which is the change in time
        :param f: array, force vector
        
        This function does not return anything
        """
        
        self.pos += np.multiply(self.vel, dt) + (np.multiply(np.divide(np.array(f), 2*self.mass), dt**2))


    def update_vel(self, dt, f):
        """
        Updates the velocity vector(vel) to the 2nd order that is within an 
        infinitesimal timestep and a force acting on the particle.
        
        :param dt: float, timestep which is the change in time
        :param f: array, force vector
        
        This function does not return anything
        """
        
        self.vel += (np.multiply(np.divide(np.array(f), 2*self.mass), dt))
                                                 

    @staticmethod
    def new_particle(file):
        """
        Initialises a Particle3D instance given an input file handle. This
        function reads in a line in the file and spilts it into list items
        where each is labels and initialised using the '__init__' method.
        
        The input file should contain one line per particle in the following format:
        <label>   <mass>  <x> <y> <z>    <vx> <vy> <vz>
        
        :param inputFile: Readable file handle in the above format

        :return Particle3D instance
        """   
        
        info = file.readline().split()
        label = str(info[0])
        mass = float(info[1])
        pos = np.array((info[2:5]), float) #Storing into np.array floats
        vel = np.array((info[5:]), float)
            
        return Particle3D(label, mass, pos, vel)


    @staticmethod
    def sys_kinetic(p3d_list):
        """
        Computes the total kinetic energy in a list of P3D's

        :param p3d_list: list in which each item is a P3D instance
        :return sys_ke: The total kinetic energy of the system
        """
        sys_ke = 0
        for i in p3d_list:
            sys_ke += i.kinetic_e()
            
        return sys_ke


    @staticmethod
    def com_velocity(p3d_list):
        """
        Computes the total mass and CoM velocity of a list of P3D's

        :param p3d_list: list in which each item is a P3D instance
        :return total_mass: The total mass of the system 
        :return com_vel: Centre-of-mass velocity
        """
        total_mass = 0
        com_vel = 0
        total_mom = 0
        
        for i in p3d_list:
            total_mass += i.mass
            total_mom += i.momentum()
    
        com_vel = total_mom / total_mass
        return total_mass, com_vel


