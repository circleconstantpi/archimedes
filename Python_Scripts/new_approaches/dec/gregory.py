#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''Calculation of the circle constant Pi.

The presented iterative algorithm bases on an geometrical approach of
James Gregory (1638 - 1675). He calculated the area of an inscribed and
and a circumscribed regular polygon to determine the value of Pi.
'''
# pylint: disable=invalid-name
# pylint: disable=line-too-long

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023, Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.1"

# Import the standard Python module math.
from decimal import Decimal as D
from decimal import getcontext, ROUND_HALF_DOWN

# Define the constants for the calculation.
PRECISION = 102
ITERATION = 166

# Set the precision and the rounding method.
getcontext().prec = PRECISION
getcontext().rounding = ROUND_HALF_DOWN

# Define the radius of the circle.
r = 1

# Define the start values.
a0 = D(6) * D(r) * (1/D(2)) * D(3).sqrt()  # outer area
b0 = D(6) * D(r) * (1/D(8)) * D(3).sqrt()  # inner area

# Run an iteration from 0 to ITERATION plus 1.
for i in range(0, ITERATION+1):
    # Use the start values in the first loop.
    if i == 0:
        a1 = a0
        b1 = b0
    else:
        # Calculate the inner and outer area.
        b1 = D(a0*b0).sqrt()
        a1 = D(2*a0*b1)/D(a0 + b1)
    # Store the old values for the next loop.
    a0 = a1
    b0 = b1
# Calculate the circle constant Pi.
# ac = (2*a1*b1) / (a1 + b1)*r  # Harmonic mean
ac = (a1 + b1) / (2*r)          # Arithmetic mean

# Print calculated and given Pi to the terminal window.
print("Calculation:", str(ac)[:102])
print("Reference:   3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679")
