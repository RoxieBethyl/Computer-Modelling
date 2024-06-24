"""
Name: Blythe Fernandes
UNN number: s1967975
"""

import numpy as np #imports numpy functions

'''
In this module, the numpy library is imported for use in the functions defined 
below.'l' defines the lenght of the box, from the origin, in which the particle
sits in. 'pos' defines the postion vector of the particle in the box.

Periodic Boundary Condiitons:
pbc(c,l) takes in the position vector(array) and is divided by the box-lenght.
The mentioned function calculates the postion of the particle if it were in the
original box of lenght l.

Minimum Image Convention:
mic(c,l) takes in the position vector(array) and calculates the position vector
of the particle that is closest to the origin. The box has a range of -l/2 to l/2
'''
def pbc(c, l):
    return np.mod(c, l)

def mic(pos, l):
    return np.mod(pos + l/2, l) - l/2

def main():
    pos = np.array([float(input("Enter the x-coodinate: ")), float(input("Enter the y-coodinate: ")), float(input("Enter the z-coodinate: "))])
    l = float(input("Enter the box lenght: "))
    while l < 0:
        print("Input postive lenght only")
        l = float(input("Enter the box lenght: "))
    
    print("\nThe position vector from the origin of the particle in the original box is: ", pbc(pos, l))
    print("\nThe position vector of the particle closest to the origin is: ", mic(pos, l))
    
if __name__ == "__main__":
    main()
