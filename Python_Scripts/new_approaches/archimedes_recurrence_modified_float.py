#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''Modern iterative Archimedes algorithm using the insights of Snellius
and the refinement from Dörrie.
'''

# Import the standard Python module math.
import math

# Define the radius.
r = 1

# Define the start values.
a0 = r * 2 * math.sqrt(3)  # half of outer perimeter
b0 = r * 3                 # half of inner perimeter

# Set the number of iterations.
iteration = 7

# Loop an iteration from 0 to 5 to get 5 values of Pi.
for i in range(0, iteration+1):
    # Use the start values in the first loop.
    if i == 0:
        a1 = a0
        b1 = b0
    else:
        # Calculate the half of inner and outer perimeter.
        a1 = (2*a0*b0)/(a0 + b0)
        b1 = math.sqrt(b0*a1)
    # Store the old values for the next loop.
    a0 = a1
    b0 = b1
    # Calculate the refinement of inner and outer bound.
    b3 = (3*a1*b1)/(2*a1 + b1)
    a3 = (a1 * b1**2)**(1/3)
    # Calculate and print the Archimedes constant.
    #ac = (a3 + b3)/2
    ac = (1/r)*(a3 + 4*b3)/5
    print(ac)

print("\n3.14159265358979323\n")
