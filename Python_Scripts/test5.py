#!/usr/bin/python3
'''Archimedes algorithm for calculating the perimeter of the
outer regular polygon based on Archimedes approach.

Description:
Pi can be calculated up to 14 correct places. Strictly spoken, this is
a upper bound of pi, since only the outer polygon is considered.

Limitation:
The code is limited to round about 1020 iterations on the test system.

Test system:
Python 3.8.10; Linux Mint 20.3 Una, Ubuntu Focal, GNU/Linux, x86_64

To-Do:
Write a version for big numbers to overcome the limitation.

Differences to SageMath:
1. The standard Python module math has to be imported.
2. Python's exponentiation operator is ** not ^.
'''
# pylint: disable=redefined-outer-name
# pylint: disable=invalid-name

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023 Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.1"

# Import the standard Python module math.
import math

# Define the function for the iterative calculation of Pi.
def archimedes_inner_polygon_method2(OA, OC, AC, iteration=5):
    '''First modified Archimedes algorithm for calculating the perimeter
    of the outer regular polygon.'''
    # Run a for loop in the range from 0 to the value of iteration.
    for i in range(0, iteration+1):
        # Calculate the number of edges.
        n = 6*2**i
        # No iteration on first loop.
        if i == 0:
            # Calculate the approximation for pi.
            ac = AC*n
        else:
            # Calculate the length of the hypotenuse and the length of the edge.
            AD = AC*OA/(OA+OC)
            OD = math.sqrt(OA**2+AD**2)
            # Store the values for the next iteration.
            AC = AD
            OC = OD
            # Calculate the approximation for pi.
            ac = AD*n
    # Return the approximation of Archimedes constant.
    return ac

# Start values for the calculation of the inner polygon.
# OB = Incircle radius
# OE = Circumcircle radius
# BE = Half of edge length
# Start values 6-gon (hexagon)
AC = 1/(math.sqrt(3))
OA = 1
OC = (2/3) * math.sqrt(3)

# Run a simple test.
iteration = 5
Pi = archimedes_inner_polygon_method2(OA, OC, AC, iteration=iteration)
print(Pi)
