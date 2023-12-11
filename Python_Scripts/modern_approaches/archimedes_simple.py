#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Import the standard Python module math.
import math

# Define the radius.
r = 1

# Define the start values.
Sn = r * (2/3) * math.sqrt(3)  # half of outer perimeter
sn = r * 1                     # half of inner perimeter

# Set the number of iterations.
iteration = 4

# Loop an iteration from 0 to 5 to get 5 values of Pi.
for i in range(0, iteration+1):
    # Calculate the number of edges.
    n = 6 * 2**i
    # Use the start values in the first loop.
    if i == 0:
        s2n = sn
        S2n = Sn
    else:
        # Calculate the half of inner and outer perimeter.
        s2n = math.sqrt(2 - math.sqrt(4 - sn**2))
        S2n = s2n / (math.sqrt(1 - (s2n/2)**2))
    # Store the old values for the next loop.
    sn = s2n
    Sn = S2n
    # Calculate and print the Archimedes constant.
    ac = (S2n + s2n)*n/4
    print(ac)

print("\n3.14159265358979323\n")
