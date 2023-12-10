#!/usr/bin/python3
'''Archimedes iterative algorithm.'''

# Import the standard Python module math.
import math

# Define the start values.
a0 = 2 * math.sqrt(3)  # half of outer perimeter
b0 = 3                 # half of inner perimeter

# Run an iteration from 0 to 5 to get 5 values of Pi.
for i in range(0, 5):
    if i == 0:
        a1 = a0
        b1 = b0
    else:
        a1 = (2*a0*b0)/(a0 + b0)
        b1 = math.sqrt(b0*a1)
    b0 = b1
    a0 = a1
    ac = (b1 + a1) / 2
    print(ac)
