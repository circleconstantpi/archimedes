#!/usr/bin/python3
'''Archimedes algorithm for calculating the perimeter of the outer regular
polygon based on Archimedes approach.

Description:
Pi can be calculated up to an infinite number of correct places. Strictly
spoken, this is a upper bound of pi, since only the inner polygon is
considered. 100 places can be reached using a precision of 204 and 167
iterations. 

Limitation:
No limitations yet.

Test system:
Python 3.8.10; Linux Mint 20.3 Una, Ubuntu Focal, GNU/Linux, x86_64

To-Do:
Nothing to do yet.

Differences to SageMath:
1. The standard Python module math has to be imported.
2. Python's exponentiation operator is ** not ^.
'''
# pylint: disable=redefined-outer-name
# pylint: disable=invalid-name

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C), 2023 Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.1"

# Import the standard Python module math.
#import math
from decimal import Decimal as D
#from decimal import getcontext, ROUND_DOWN
from decimal import *

# Set precision.
getcontext().prec = 204
getcontext().rounding = ROUND_DOWN

# Define the function for the iterative calculation of Pi.
def outer_from_inner_polygon(AB, AC, BC, iteration=5, verbose=True):
    '''Archimedes algorithm for calculating the perimeter
    of the outer regular polygon.'''
    # Run a for loop in the range from 0 to the value of iteration.
    for i in range(0, iteration+1):
        # Calculate the number of edges.
        n = 6*2**i
        # No iteration on first loop.
        if i == 0:
            # Calculate the approximation for pi.
            ac = (D(BC)/D(AC))*n
        else:
            # Calculate the length of the hypotenuse and the length of the edge.
            AD = D(AB)/D((BC**2/(AB + AC)**2) + 1).sqrt()
            BD = D(AB**2 - AD**2).sqrt()
            # Store the values for the next iteration.
            BC = D(BD)
            AC = D(AD)
            # Calculate the approximation for pi.
            ac = (D(BD)/D(AD))*n
    # Return the approximation of Archimedes constant.
    return ac

# Start values 6-gon (hexagon) for the calculation of the inner polygon.
# OB = Incircle radius
# OE = Circumcircle radius
# BE = Half of edge length
AC = D(3).sqrt()
AB = D(2)
BC = D(1)

# Run a simple test.
iteration = 167
Pi = outer_from_inner_polygon(AB, AC, BC, iteration=iteration)
print("{0:.110f}".format(Pi))
print("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679")
