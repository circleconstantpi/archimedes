#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''Archimedes algorithm.

The programmed iterative Archimedes algorithm uses the approach
from Dörrie and the refinement from Netz for the improvement of
the calculation result of Pi.
'''
# pylint: disable=invalid-name

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023, Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.1"

# Import the standard Python module math.
from decimal import Decimal as D
from decimal import getcontext, ROUND_HALF_DOWN

# Initialise the constants.
PRECISION = 102
ITERATION = 54

# Set the precision and the rounding method.
getcontext().prec = PRECISION
getcontext().rounding = ROUND_HALF_DOWN

# Define a heredoc consisting of Pi with 100 places.
PI100 = '''
3.
14159265358979323846264338327950288419716939937510
58209749445923078164062862089986280348253421170679
'''

# ----------------------------------------------------------------------
# Helper function remove_ws()
# ----------------------------------------------------------------------
def remove_ws(instr):
    '''Remove whitespaces defined by a list.
    '''
    # Define the whitespaces to remove.
    mapping = [("\n", ""), ("\r", ""), ("\t", ""), (" ", "")]
    # Remove the whitespaces.
    for k, v in mapping:
        instr = instr.replace(k, v)
    # Return trimmed string.
    return instr

# ----------------------------------------------------------------------
# Helper function correct_digits()
# ----------------------------------------------------------------------
def correct_digits(chkpi, refpi=''):
    '''Calculate the correct digits of a given pi number.'''
    # Define three parts of circle number pi.
    a = "3."
    b = "14159265358979323846264338327950288419716939937510"
    c = "58209749445923078164062862089986280348253421170679"
    if refpi == '':
        refpi = a + b + c
    if len(str(refpi)) < len(str(chkpi)):
        return None, None
    correct = ''
    idx = 0
    for char in str(chkpi):
        if char == str(refpi)[idx]:
            correct += char
            idx += 1
        else:
            break
    return (correct, idx-2)

# Create PI100.
PI100 = remove_ws(PI100)

# Define the radius.
r = 1

# Define the start values.
a0 = D(r) * D(2) * D(3).sqrt()   # half of outer perimeter
b0 = D(r) * D(3)                 # half of inner perimeter

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
    # Calculate the refinement of inner and outer bound.
    b3 = D(3*a1*b1)/D(2*a1 + b1)
    a3 = D(a1 * b1**2)**(D(1)/D(3))
# Calculate and print the Archimedes constant.
ac = D((D(1)/D(r))*D(a3 + 4*b3))/D(5)

print("Calculation:", str(ac)[:102])
print("Reference:  ", PI100)

circle_constant, number_digits = correct_digits(str(ac)[:102], PI100)
print("\nMatching places:", number_digits, " -> ", circle_constant)
