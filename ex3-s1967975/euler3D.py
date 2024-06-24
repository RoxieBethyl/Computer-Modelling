"""
Author: Blythe Fernandes
UUN number: s1967975

 CompMod Ex3: symplectic Euler time integration of a particle moving with the 
 morse potential.

 Produces plots of the difference between the positions of the particles and 
 its energy, both as function of time. Also saves both to file.

 The potential is U_M(x) = D_e * {[1 - exp(-α(r_12-r_e))] ^ 2 - 1}, where D_e, 
 α, r_e are read in the main() method and passed to the functions that 
 calculates force and potential energy.
"""

import sys
import numpy as np
import matplotlib.pyplot as pyplot
from particle3D import Particle3D as p3d
from scipy.signal import find_peaks
from scipy.constants import c

def force_dw(D, r_e, alpha, r_12):
    """
    Method to return the force on particle 1 that expriences a morse potential.
    Force is given by
    F(x) = 2*α*D_e * { [1-exp(-α(r_12-r_e))] * exp(-α(r_12-r_e)) } * r̂_12

    :param D: Constant of the system
    :param r_e: Constant of the system
    :param alpha: Constant of the system
    :param r_12: Distance between to particles over time
    :return: force acting on particle as Numpy array
    """
    
    r_12n = np.linalg.norm(r_12)
    r_12hat = r_12/r_12n
    return (2 * alpha * D * (1-np.exp(-alpha*(r_12n-r_e))) * np.exp(-alpha*(r_12n-r_e))) * r_12hat


def pot_energy_dw(D, r_e, alpha, r_12):
    """
    Method to return the potential energy of 2 particles of a morse potential.
    Potential is defined as 
    U_M(x) = D_e * {[1 - exp(-α(r_12-r_e))] ^ 2 - 1}

    :param D: Constant of the system
    :param r_e: Constant of the system
    :param alpha: Constant of the system
    :param r_12: Distance between to particles over time
    :return: force acting on particle as Numpy array
    """
    
    return (D * ((1 - (np.exp( -alpha * (np.linalg.norm(r_12)-r_e))))**2 - 1))


def frequency(x, t):
    """
    Method to return the frequency of the funtion and returns this value.
    Uses scipy funtion to find the peaks of the data. Converts the period in 
    the unit used in the simulation to calculate the frequency.
    
    :param x: Particle position difference (array)
    :param t: Time (array)
    
    :retrun f: The frequency of the system
    """
    
    peak = find_peaks(x)[0]
    period = (t[peak[2]] - t[peak[1]]) * 1.018050571e-14

    return 1/period


# Begin main code
def main():
    """
    A simulation is conducted by using the symplectic Euler algorithm which 
    reads input from a given file and outputs the data into another file.
    This simulation updates the using the 2nd differential equtions of the 
    position.
    
    This function prints 2 graphs that show the energy vs time and the 
    difference between the 2 particle in this simulation against time.
    
    A while-loop construct is used to read input data from the file and set the
    initial parameter. These are inialised for the Particle3D (p3d) class to 
    make calculations and gather data.
    
    This data is used to find:
     - The distance between the parameters
     - The total energy between the particle (using the sys_kinetic, from p3d 
                                              and pot_energy_dw methods)
    
    A For-loop is used to perform the time integration where the postions, 
    velocities, force and the energies are updated. The frequency for the 
    simulation is calculated.
    """
    
    # Read name of output file from command line
    if len(sys.argv)!=2:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <output file>")
        sys.exit()
    else:
        outfile_name = sys.argv[1]

    outfile = open(outfile_name, "w")     # Open output file 
    
    
    
    # Conditional statement to write out initial conditions from file
    time = 0          # Set up simulation parameter for time
    with open(input("Enter file name: ")) as f:
        # Set up particle initial conditions (intialised in this order):
        #  label - string
        #  mass - float
        #  pos - np.array[type=float]
        #  vel - np.array[type=float]
        
        f.readline()    # skips comment line
        
        #Data from line 1
        line1 = f.readline().split()
        dt = float(line1[0])
        numstep = int(line1[1])
        
        f.readline()    # skips comment line
        
        #Data from line 2
        line2 = f.readline().split()
        D = float(line2[0])
        r_e = float(line2[1])
        alpha = float(line2[2])
        
        #Data from line 3
        line3 = f.readline().split()
        p1 = p3d(str(line3[0]), float(line3[1]), np.array(line3[2:5],float), np.array(line3[5:], float))
        
        #Data from line 4
        line4 = f.readline().split()
        p2 = p3d(str(line4[0]), float(line4[1]), np.array(line4[2:5],float), np.array(line4[5:], float))
        
        
        
    # Write out initial conditions
    r_12 = np.subtract(p2.pos, p1.pos)      # Distance between the particals
    particle_list = [p1, p2]                # Writes list
    i_energy = p3d.sys_kinetic(particle_list) + pot_energy_dw(D, r_e, alpha, r_12)
                                            # Get initial energy



    # Writing values(time, particle seperation and energy) to file
    outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time, np.linalg.norm(r_12), i_energy))
    
    
    
    # Initialise data lists for plotting later
    time_list = [time]
    diff_list = [np.linalg.norm(r_12)]
    energy_list = [i_energy]



    # Start the time integration loop
    for i in range(numstep):
        r_12 = np.subtract(p2.update_pos(dt), p1.update_pos(dt))    # Update particle position
        force = force_dw(D, r_e, alpha, r_12)       # Calculate force
        
        # Update particle velocity 
        p1.update_vel(dt, force)
        p2.update_vel(dt, -force)   
        
        particle_list = [p1, p2]    # Rewrites list with updates
        
        time += dt  # Increase time
        
        # Calculates energy for this update instance
        energy = p3d.sys_kinetic(particle_list) + pot_energy_dw(D, r_e, alpha, r_12)
        # Output particle information
        outfile.write("{0:f} {1:f} {2:12.8f}\n".format(time, np.linalg.norm(r_12), i_energy))

        # Append information to data lists
        time_list.append(time)
        diff_list.append(np.linalg.norm(r_12))
        energy_list.append(energy)
    
        if energy > i_energy: 
            E_max = float(energy)
            
        if energy < i_energy: 
            E_min = float(energy)
            
        
    
    # Post-simulation:
    outfile.close()  # Close output file
    
    # Plot particle difference over time to screen
    pyplot.title('Symplectic Euler: Difference vs Time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Particle position difference')
    pyplot.plot(time_list, diff_list)
    pyplot.show()

    # Plot particle energy to screen
    pyplot.title('Symplectic Euler: Total Energy vs Time')
    pyplot.xlabel('Time')
    pyplot.ylabel('Energy')
    pyplot.plot(time_list, energy_list)
    pyplot.show()

    f = frequency(diff_list, time_list)
    
    # Prints results
    print("Energy inaccuracy: ", ((E_max - E_min)/(i_energy)))
    print("Frequency: ", f)
    print("Wave number: ", (f/(c*(100)))) # Diving f by c



# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
