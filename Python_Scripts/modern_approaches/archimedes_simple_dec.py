#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''Modern iterative Archimedes algorithm using the insights of Snellius
and the refinement from Dörrie.
'''

# Import the standard Python module math.
from decimal import Decimal as D
from decimal import getcontext

# Set precision.
getcontext().prec = 253

# Define the radius.
r = 1

# Define the start values.
Sn = r * (D(2)/D(3)) * D(3).sqrt()  # outer edge
sn = r * 1                          # inner edge

# Set the number of iterations.
iteration = 252

# Loop an iteration from 0 to 5 to get 5 values of Pi.
for i in range(0, iteration+1):
    # Calculate the number of edges.
    n = 6 * 2**i
    # Use the start values in the first loop.
    if i == 0:
        s2n = D(sn)
        S2n = D(Sn)
    else:
        # Calculate the half of inner and outer perimeter.
        s2n = D(2 - D(4 - sn**2).sqrt()).sqrt()
        S2n = s2n / (1 - (s2n/2)**2).sqrt()
    # Store the old values for the next loop.
    sn = D(s2n)
    Sn = D(S2n)
# Calculate and print the Archimedes constant.
ac = D(S2n + s2n)*n/D(4)
print("Calculation: {}".format((str(ac))[:102]))

print("Reference:   3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679")
