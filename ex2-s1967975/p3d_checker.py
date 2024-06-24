#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file imports and checks the Particle3D class by performing some basic
tests. Amongst the things tested are:
    * basic instantiation
    * instance methods
    * creating instances from a data file
    * static methods

@author: Blythe Fernandes
"""
import sys
import numpy as np
import math
from particle3D import Particle3D as p3d

tests_failed = 0

# Basic instance, using arrays
m0 = 1.0
x0 = np.zeros(3)
v0 = np.array([1, -1, 0], float)
try:
    new_p3d = p3d('new', m0, x0, v0)
except:
    print("Cannot use __init__. Aborting...")
    sys.exit(-1)

# Test 1: kinetic energy
ke0 = 1.0
ke = new_p3d.kinetic_e()
if math.isclose(ke0, ke):
    print("Kinetic energy check 1: ok")
else:
    print("Kinetic energy check 1: ERROR")
    tests_failed += 1

# Test 3: Updaters. Resets pos/vel in between runs.
dt = 0.125
force = np.ones(3)
dx1 = 0.125 * np.array([1, -1, 0], float)

# up1
new_p3d.update_pos(dt)
if np.allclose(dx1, new_p3d.pos):
    print("Position updater 1 check: ok")
else:
    print("Position updater 1 check 1: ERROR")
    tests_failed += 1
new_p3d.pos = x0


# Test the file handle creator
with open('data_for_p3d') as f:
    f.readline()  # prune comment line
    alice = p3d.new_particle(f)
    bob = p3d.new_particle(f)
    p3d_list = [alice, bob]

# And now we have a list of p3d's, check a static method...
total_ke0 = 0.625
total_ke = p3d.sys_kinetic(p3d_list)
if math.isclose(total_ke0, total_ke):
    print("Kinetic energy check 2: ok")
else:
    print("Kinetic energy check 2: ERROR")
    tests_failed += 1
    
# Check __str__
print("\nChecking __str__. Your particle output is:\n{0:s}".format(str(alice)))
print("\nIt should be similar to:\nAlice    0.0   0.0   0.0\n")

print(f"Tests completed. There were {tests_failed} errors.")

