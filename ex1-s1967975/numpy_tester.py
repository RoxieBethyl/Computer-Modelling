"""
Name: Blythe Fernandes
UNN number: s1967975
"""

'''
This modules immports the numpy and random library to create random vectors
to test the mentioned vector identities using numpy functions defined in the
numpyFunctions.py module
'''

import random as rnd
import numpyFunctions as npy
import numpy as np

# Main method
def main():

    #Test and check indentites using funtions from npyectors.py
    v1 = np.array([rnd.random(), rnd.random(), rnd.random()])
    v2 = np.array([rnd.random(), rnd.random(), rnd.random()])
    v3 = np.array([rnd.random(), rnd.random(), rnd.random()])
    
    #Print out all npyectors
    print("v1 = ", v1)
    print("v2 = ", v2)
    print("v3 = ", v3)
    
    #Checks [v1 x v2 = -v2 x v1] indentity
    print("\n")
    npy.same(npy.cross(v1,v2), npy.cross(npy.multi_scale(v2,-1),v1))
    
    #Checks [v1 x (v2 + v3) = (v1 x v2) + (v1 x v3)]
    print("\n")
    npy.same(npy.cross(v1, npy.add(v2, v3)), npy.add(npy.cross(v1, v2),npy.cross(v2, v3)))
    
    #Checks [v1 x (v2 x v3) = (v1 . v3)*v2 - (v1 . v2)*v3]
    print("\n")
    npy.same(npy.cross(v1, npy.cross(v2, v3)),npy.diff(npy.multi_scale(v2, npy.dot(v1, v3)),npy.multi_scale(v3, npy.dot(v1, v2))))
        
main()