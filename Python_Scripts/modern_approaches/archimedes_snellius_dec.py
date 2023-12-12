#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''Modern iterative Archimedes algorithm using the insights with
respect too the convergence of the inner and the outer bound from
Snellius and the refinement of the lower and the upper bound from
Dörrie.
'''

# Import the standard Python module math.
from decimal import Decimal as D
from decimal import getcontext, ROUND_HALF_EVEN

# Define the constants for the calculation.
PRECISION = 102
ITERATION = 83

# Set the precision and the rounding method.
getcontext().prec=PRECISION
getcontext().rounding=ROUND_HALF_EVEN

# Define the radius.
r = 1

# Define the start values.
a0 = r * 2 * D(3).sqrt()   # half of outer perimeter
b0 = r * 3                 # half of inner perimeter

# Set the number of iterations.

# Loop an iteration from 0 to ITERATION plus 1.
for i in range(0, ITERATION+1):
    # Use the start values in the first loop.
    if i == 0:
        a1 = D(a0)
        b1 = D(b0)
    else:
        # Calculate the half of inner and outer perimeter.
        a1 = D(2*a0*b0)/D(a0 + b0)
        b1 = D(b0*a1).sqrt()
    # Store the old values for the next loop.
    a0 = D(a1)
    b0 = D(b1)
# Calculate and print the Archimedes constant.
ac = D((D(1)/D(r))*D(a1 + 2*b1))/D(3)

# Print calculated and given Pi to the terminal window.
print("Calculation:", str(ac)[:102])
print("Reference:   3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679")
