"""
Name: Blythe Fernandes
UNN number: s1967975
"""

'''
This module imports the numpy and random library where it's functions are used
to test the vector identities.
'''

import random as rnd
import vectors as V

# Main method
def main(): #Test and check indentites using funtions from vectors.py

    #Sets vectors used to test the identites in this module
    v1 = [rnd.random(), rnd.random(), rnd.random()]
    v2 = [rnd.random(), rnd.random(), rnd.random()]
    v3 = [rnd.random(), rnd.random(), rnd.random()]
    
    #Print out all vectors
    print("v1 = ", v1)
    print("v2 = ", v2)
    print("v3 = ", v3)
    
    
    #Checks [v1 x v2 = -v2 x v1] indentity
    print("\n")
    V.same(V.cross(v1,v2), V.cross(V.multi_scale(v2,-1),v1))
    
    #Checks [v1 x (v2 + v3) = (v1 x v2) + (v1 x v3)]
    print("\n")
    V.same(V.cross(v1, V.add(v2, v3)), V.add(V.cross(v1, v2),V.cross(v2, v3)))
    
    #Checks [v1 x (v2 x v3) = (v1 . v3)*v2 - (v1 . v2)*v3]   
    print("\n")
    V.same(V.cross(v1, V.cross(v2, v3)),V.diff(V.multi_scale(v2, V.dot(v1, v3)),V.multi_scale(v3, V.dot(v1, v2))))
    
main()