"""
Author: Blythe Fernandes
UUN Number: s1967975

 Project B: Lennard-Jones system
 This is the main code file to describe an N-body system interacting with the
 the Lennard-Jones pair potential. The code shall, eventually, be able to
 simulate the solid, liquid, and gas phases of Argon using periodic boundary
 conditions(pbc) and the minimum image covention(mic). From these the 
 equilibrium properties can be found for all these states. Functions such as
 are used to shorten main code.
 
"""

import numpy as np
from numpy import linalg as lin
import mdutilities as mdu
from particle3D import Particle3D as p3d
import matplotlib.pyplot as pyplot
import pbc
import copy


def calc_diff(particles, n, boxSize):
    """
    Function to calculate all the seperations for all particles with all the
    other particles.

    Parameters:
        particles(list): List of all particle data initialised using the
                         Particle3D class
        n(integer): Number of particles in particle list
        boxSize(float): Length of the box
        
    Returns:
        diff_list(np.array): List of all the seperations between all the
                             particles
    """

    diff_list = np.zeros((n, n, 3))
    for i in range(n):
        for j in range(n):
              if i > j:
                  diff_list[i, j] = pbc.mic(particles[i].position -
                                            particles[j].position, boxSize)
                  diff_list[j, i] = -diff_list[i, j]

    return diff_list


def calc_potential(n, diff_list, a):
    """
    Function to calculate the total potential energy of each particle
    interaction.

    Parameters:
        n (int): Number of particles
        diff_lisy (3d np.array): The difference between the particles within 
                                 the system
                                 
    Returns:
        pot(float): The total sum of tpotential energy of the system.
    """
    pot = 0
    for i in range(n):
        for j in range(n):
            r = lin.norm(diff_list[i, j])
            if i > j and r <= a:                
                r = lin.norm(diff_list[i, j])
                pot += 4. * ((r**-12) - (r**-6))
                
    return pot


def calc_force(n, diff_list, a):
    """
    Function to calculate the total force beteween one particle and the rest
    of the particles in the system.

    Parameters:
        n (int): Number of particles
        diff_lisy (3d np.array): The difference between the particles within 
                                 the system

    Returns:
        force_list(np.array): The total sum of the forces on the particle
    """
    force_list = np.zeros((n, n, 3))
    for i in range(n):
        for j in range(n):
            r = lin.norm(diff_list[i, j])
            if i > j and r <= a:
                force_list[i, j] = diff_list[i, j] * -48. * ((r**-14) - (0.5 * 
                                                                     (r**-8)))
                force_list[j, i] = -force_list[i, j]
                
    return np.sum(-force_list, axis = 1)


def update_pos(particles, n, dt, force_list, boxSize):
    """
    Loop to calculate update all the position simultaneouly between each
    particles.

    Parameters:
        particle (3darray): Array of particle data
        n (int): Number of particles
        dt(float): time step
        force_list(np.array): The total force of the system
        boxSize(float): Length of the box
    """
    
    for i in range(n):
        p3d.update_pos_2nd(particles[i], dt, force_list[i])
        particles[i].position = pbc.pbc(particles[i].position, boxSize)

def update_vel(particles, n, dt, force, force_new):
    """
    Loop to calculate update all the velocties due to each particle
    interaction.

    Parameters:
        particle (3darray): Array of particle data
        n (int): Number of particles
        dt(float): time step
        force(np.array): The total force of the system
        force_new(np.array): The new total force of the system
    """
    for i in range(n):
        p3d.update_vel(particles[i], dt, (0.5 * (force[i] + force_new[i])))
        

def rdf(diff_list, n, numstep, rho):
    """
    Function to find the probability to find a particle at the various 
    distances.

    Parameters:
        diff_data(3d np.array): Data of each particle seperation
        n (int): Number of particles
        numstep(int): The time step at which the particle is at
        rho(float): The density of the material
        
    Returns:
        distance(np.array): Array storing data of average particle movements
        g(np.array): Array of a the RDFs found at that timestep
    """
    
    r = lin.norm(diff_list, axis = -1)
    index, blocks = np.histogram(r, bins = 100, range = (0,5))
    distance = blocks[1:-1]
    g = index[1:]/(distance + 0.5*(blocks[1] - blocks[0]))**2
    
    return distance, g

def main():
    """
    Particles are choosen to completely fill a lattice the number of particles
    should be 4n^3 for some integer n. In this loop, the particle positions 
    are updated and the particle seperations between all the particles in the 
    system are stored. 
    
    Parameters are initialised using functions and the time integration loop is
    started. The particle position, velocities are initialised using the
    provided mdutilities module. 
    
    The MSD, RDF, forces and energies are using to analyse the Lennard-Jones
    System.
    """
    
    particles = []
    
    '''
    # Test 1
    r_0 = np.power(2., 1./6.)
    n = 2
    dt = 0.01
    numstep = 1000
    boxSize = np.array([3., 3., 3.])
    
    # initialise particles with particle3D as 'label, mass, pos(xyz), vel(xyz)'
    for i in range(1, n+1):
        # For odd particles
        if i%2 != 0:
            # input("Enter the position vector: ")
            particles += [p3d(str(i), 14., pbc.pbc(np.array([0., 0., 0.]), 
                                          boxSize), np.array([0.01, 0., 0.]))]

        # Even
        else:
            particles += [p3d(str(i), 14., pbc.pbc(np.array(
                   [r_0 + 0.1, 0., 0.]), boxSize), np.array([-0.01, 0., 0.]))]
    '''
    
    # Test 2q
    # User input for the number of particles, density, temperature and timestep
    # n = 108           #int(input("Enter the number of particles: "))
    # den = 0.0843         #float(input("Enter the density: "))
    # T = 0.6          #float(input("Enter the temperature: "))
    # dt = 0.001        #float(input("Enter time increment: "))
    # numstep = 10000   #int(input("Enter the timestep: "))
    # a = 3.5          #float(input("Enter the cut-off radius: "))
    
    n = 30
    den = 0.1
    T = 1.
    dt = 0.01
    numstep = 1000
    a = 3.5
    
    # Create particles
    for i in range(n):
        # Initialise particles with particle3D as 'label, mass, pos(xyz),
        # vel(xyz)'
        particles += [p3d(str(i), 1., np.array([0., 0., 0.]), 
                          np.array([0., 0., 0.]))]

    # Intialise particles for position and velocity
    boxSize, lattice = mdu.set_initial_positions(den, particles)  
    mdu.set_initial_velocities(T, particles)

    # Setting intial paramenters
    time = 0
    r = 0.
    g = 0.
    time_list = []
    energy_list = []
    kinetic = []
    potential = []
    g_data = 0
    diff_data = np.zeros((numstep, n, n, 3))
    msd_data = np.zeros((numstep + 1))
    
    # Initalising seperations between all particles with each other
    force_list = calc_force(n, calc_diff(particles, n, boxSize[0]), a)
    kinetic.append(p3d.sys_kinetic(particles))
    potential.append(calc_potential(n, calc_diff(particles, n, boxSize[0]), a))
    energy_list.append(kinetic[0] + potential[0])
    time_list.append(time)
    
    # Making a copy of the initial position for msd
    i_particles = copy.deepcopy(particles)
    
    # Calculating the Mean Square Displacement at the zeroth time step 
    for i in range(n):
        msd_data[0] += (np.linalg.norm(pbc.mic((particles[i].position - 
                                        i_particles[i].position), boxSize))**2)/n
                
    # Integration Loop order:
    #   Update particle position
    #   Seperations found
    #   LJ force force
    #   Update particle velocites
    #   Rewrite force
    #   Energies found
        
    # Start the time integration loop
    # Opening files to write Energies, MSDs and RDFs to file
    with open('outfile.xyz', 'w') as outfile, open('msd.dat', 'w') as msdFile,\
    open('rdf.dat', 'w') as rdfFile:
        msdFile.write("Time MSD\n")
        for t in range(numstep):
            
            # Output particle information as XYZ file to use in VMD
            outfile.write("{0}\n" .format(n))
            outfile.write("Point = {0}\n" .format(t))
            for i in range(len(particles)):
                outfile.write("{0} {1} {2} {3}\n" .format(i, particles[i]
                .position[0], particles[i].position[1], particles[i]
                .position[2]))
            
            # Calculating the MSD
            for i in range(n):
                msd_data[t+1] += (np.linalg.norm(pbc.mic((particles[i].position - 
                                        i_particles[i].position), boxSize))**2)/n
            
            # Writing MSD information to file
            msdFile.write("{0} {1}\n" .format(t, msd_data[t+1]))
            
            update_pos(particles, n, dt, force_list, boxSize[0])
            diff_data[t] = calc_diff(particles, n, boxSize[0])
            
            # Update force
            new_force_list = calc_force(n, diff_data[t], a)
    
            # Update particle velocity by averaging current and new forces
            update_vel(particles, n, dt, force_list, new_force_list)
            
            # Re-define force values
            force_list = new_force_list
    
            # Calculates energy for this update instance
            kinetic.append(p3d.sys_kinetic(particles))
            potential.append(calc_potential(n, diff_data[t], a))
            energy_list.append(kinetic[t] + potential[t])
            
            r, g = rdf(diff_data[t], n, numstep, den)
            g_data += g
            # Output RDF information to file
            rdfFile.write("{0}\n {1}\n {2}\n\n" .format(t, r, g_data))
               
            # Increase time
            time += dt
    
            # Append information to data lists
            time_list.append(time)
            
            if ((t/numstep)*100) % 10 == 0:
                print((t/numstep)*100,"%")
            
    g_data = g_data/(numstep*n*4*np.pi*r**2*den)
    meanEtot = np.mean(energy_list, axis = 0)
    
    print("Mean Energy: ", meanEtot)
    
    # Post-simulation:
    # Closing all open files
    outfile.close()
    msdFile.close()
    rdfFile.close()
    
    # Plotting Graphs for the total energy, MSD and RDF
    pyplot.title('Kinetic Energy')
    pyplot.xlabel('Time(√(m/Ɛ))')
    pyplot.ylabel('Kinetic Energy(Ɛ)')
    pyplot.plot(time_list, kinetic, color = "brown")
    pyplot.show()
    
    pyplot.title('Potential Energy')
    pyplot.xlabel('Time(√(m/Ɛ))')
    pyplot.ylabel('Potential Energy(Ɛ)')
    pyplot.plot(time_list, potential, color = "green")
    pyplot.show()
    
    pyplot.title('Total Energy')
    pyplot.xlabel('Time(√(m/Ɛ))')
    pyplot.ylabel('Total Energy(Ɛ)')
    pyplot.plot(time_list, energy_list, color = "teal")
    pyplot.show()
    
    pyplot.title('Mean Square Displacement')
    pyplot.xlabel('Time(√(m/Ɛ))')
    pyplot.ylabel('MSD(σ²)')
    pyplot.plot(time_list, msd_data, color = "red")
    pyplot.show()
    
    pyplot.title('Radial Distribution Function')
    pyplot.xlabel('Distance(r)')
    pyplot.ylabel('RDF')
    pyplot.plot(r, g_data, color = "orange")
    pyplot.show()
    
# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()