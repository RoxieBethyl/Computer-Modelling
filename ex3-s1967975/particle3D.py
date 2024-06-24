"""
Author: Blythe Fernandes
UUN number: s1967975

 CompMod Ex3: Particle3D, a class to describe point particles in 3D space

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
        self.pos = np.array(pos, float) #Changing position vector into array float using numpy array
        self.vel = np.array(vel, float) #Changing velocity vector into array float using numpy array
        

    def __str__(self):
        """        
        Define output format.
        :return (string) format: <label>, pos = <[x y z]>, v = <[v_x v_y v_z]>, m = <mass>
        For particle Oxygen (2.0, 0.5, 0.3) (2.0, 0.5, 0.3) 16.0 this will print as
        "Oxygen, pos = [2.0 0.5 0.3],  v = [2.0 0.5 0.3], m = 16.0"
        """
        return str(self.label) + ", pos = " + str(self.pos) + ", v = " + str(self.vel) + ", m = " + str(self.mass)

    
    def kinetic(self):
        """
        Returns the kinetic energy of a Particle3D instance

        :return ke: float, 1/2  * m * v**2
        """
        return (0.5 * self.mass * np.linalg.norm(self.vel)**2)


    def momentum(self):
        """
        Returns the momentum of the Particles
        
        :return p: vector, m * vel(vector)
        """
        #Using a numpy function
        return np.multiply(self.vel, self.mass)


    def update_pos(self, dt):
        """
        Updates the position vector(pos) to the 1st order that is within an 
        infinitesimal timestep for the change in the particles position.
        
        :param dt: float, timestep which is the change in time
        
        :return pos: Updated position for Particle3D instance
        """
        self.pos += np.multiply(self.vel, float(dt))
        return self.pos


    def update_pos_2nd(self, dt, f):
        """
        Updates the position vector(pos) to the 2nd order that is within an 
        infinitesimal timestep and a force acting on the particle.
        
        :param dt: float, timestep which is the change in time
        :param f: array, force vector
        
        :return pos: Updated position for Particle3D instance
        """
        self.pos += np.multiply(self.vel, float(dt)) + np.multiply(np.divide(f, self.mass), 0.5 * float(dt)**2)
        return self.pos


    def update_vel(self, dt, f):
        """
        Updates the velocity vector(vel) to the 2nd order that is within an 
        infinitesimal timestep and a force acting on the particle.
        
        :param dt: float, timestep which is the change in time
        :param f: array, force vector
        
        :return vel: Updated velocity for Particle3D instance
        """
        
        self.vel += np.multiply(np.divide(f, self.mass), float(dt))
        return self.vel                                        


    @staticmethod
    def new_particle(file):
        """
        Initialises a Particle3D instance given an input file handle. This
        function reads in a line in the file and spilts it into list items
        where each is labels and initialised using the '__init__' method.
        
        The input file should contain one line per particle in the following format:
        <label>   <mass>  <pos> <y> <z>    <vx> <vy> <vz>
        
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
            sys_ke += i.kinetic()
            
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


