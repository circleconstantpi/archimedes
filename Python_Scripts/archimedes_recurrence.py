#!/usr/bin/python3

# Import the standard Python module math.
import math

# Define the startvalues.
a0 = 3
b0 = 2*math.sqrt(3)

# Run an iteration from 0 to 5 to get 5 values of Pi.
for i in range(0, 5):
    if i == 0:
        b1 = b0
        a1 = a0
    else:
        b1 = (2*b0*a0) / (b0 + a0)
        a1 = math.sqrt(b1 * a0)
    b0 = b1
    a0 = a1
    ac = (b1 + a1) / 2
    print(ac)
