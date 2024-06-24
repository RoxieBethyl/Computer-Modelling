"""
Name: Blythe Fernandes
UNN number: s1967975
"""

import math

def sq_mod(v):
    '''
    Square Modulus of vector (3 components)
    
    :param v: vector (v1, v2, v3)
    :return: v1^2 + v2^2 + v3^2
    '''
    return v[0]**2 + v[1]**2 + v[2]**2


def sqrt_mod(v):
    '''
    Square root Modulus vector
    
    :param v: vector (v1, v2, v3)
    :Return: sqrt(v1^2 + v2^2 + v3^2)
    '''
    return math.sqrt(sq_mod(v))

def multi_scale(v, scalar):
    '''
    Multiplication of vector with a scalar
    
    :param v: vector (v1, v2, v3)
    :param scalar: scaling factor
    :return: scaled vector (v1*scalar, v2*scalar, v3*scalar)
    '''
    return [v[0]*scalar, v[1]*scalar, v[2]*scalar]

def div_scale(v, factor):
    '''
    Division of vector by a factor
    
    :param v: vector (v1, v2, v3)
    :param factor: reducing factor
    :return: scaled vector (v1/factor, v2/factor, v3/factor)
    '''
    return [v[0]/factor, v[1]/factor, v[2]/factor]

def add(v1, v2):
    '''
    Vector Addition
    
    :param v1: First Vector
    :param v1: Second Vector
    :return: v1(a1, b1, c1) + v2(a2, b2, c2) = vSum(a1+a2, b1+b2, c1+c2)
    '''
    return [v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]]

def diff(v1, v2):
    '''
    Vector Difference
    
    :param v1: First Vector
    :param v1: Second Vector
    :return: v1(a1, b1, c1) - v2(a2, b2, c2) = vSub(a1-a2, b1-b2, c1-c2)
    '''
    return [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]]

def cross(v1, v2):
    '''
    Cross Product of Vectors
    
    :param v1: First Vector
    :param v2: Second Vector
    :return: the determint of the vectors 
             [i((a2*b3)-(a3*b2)), j((a3*b1)-(a1*b3)), k((a1*b2)-(a2*b1))]
    '''
    return [((v1[1]*v2[2])-(v1[2]*v2[1])), ((v1[2]*v2[0])-(v1[0]*v2[2])), ((v1[0]*v2[1])-(v1[1]*v2[0]))]

def dot(v1, v2):
    '''
    Dot Product of Vectors
    
    :param v1: First vector
    :param v2: Second Vector
    :return: the dot product, (a1*b1) + (a2*b2) + (a3*b3)
    '''
    return ((v1[0]*v2[0]) + (v1[1]*v2[1]) + (v1[2]*v2[2]))

def same(v1,v2):
    '''
    Checks for the equality of the indentity
    
    :param v1: First vector
    :param v2: Second Vector
    :return: The vectors are the same if true. If not true then then aren't the same.
    '''
    print("The 2 vectors to check are: ", v1, " and ", v2)
    
    if math.isclose(v1[0], v2[0], abs_tol=1e-09) and math.isclose(v1[1], v2[1], abs_tol=1e-09) and math.isclose(v1[2], v2[2], abs_tol=1e-09):
        print("These vectors are the same")
    else:
        print("These vectors are not the same")
        
    return