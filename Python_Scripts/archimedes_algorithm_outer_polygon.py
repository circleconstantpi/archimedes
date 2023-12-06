#!/usr/bin/python3
'''Archimedes algorithm for calculating the perimeter of the
outer regular polygon based on Archimedes approach.

Description:
Pi can be calculated up to 14 correct places. Strictly spoken, this
is a upper bound of Pi, since only the outer polygon is considered.

Limitation:
The code is limited to round about 1020 iterations on the test system.

Test system:
Python 3.8.10; Linux Mint 20.3 Una, Ubuntu Focal, GNU/Linux, x86_64

To-Do:
Write a version for big numbers to overcome the limitation.

Differences to SageMath:
1. The standard Python module math has to be imported.
2. Python's exponentiation operator is ** not ^.

See also:
docs.python.org/3.8/tutorial/floatingpoint.html
docs.python.org/3.8/library/fractions.html
docs.python.org/3.8/library/decimal.html
'''
# pylint: disable=redefined-outer-name
# pylint: disable=invalid-name

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023, Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.3"

# Import the standard Python module math.
import math

# Define the function for the iterative calculation of Pi.
def archimedes_outer_polygon(OA, OC, AC, iteration=5):
    '''Archimedes algorithm for calculating the perimeter of the outer
    regular polygon.'''
    # Store incircle radius for later use.
    r = OA
    # Run a for loop in the range from 0 to the value of iteration plus 1.
    for i in range(0, iteration+1):
        # Calculate the number of edges.
        n = 6*2**i
        # No iteration on first loop.
        if i == 0:
            # Calculate the approximation for pi.
            ac = (AC/r)*n
        else:
            # Catch an overflow error. Returns the last valid value.
            try:
                # Calculate the length of the hypotenuse and the length of the edge.
                AD = AC*OA/(OA+OC)
                OD = math.sqrt(OA**2+AD**2)
                # Store the values for the next iteration.
                AC = AD
                OC = OD
                # Calculate the approximation for pi.
                ac = (AD/r)*n
            except OverflowError:
                # Print an error message.
                print("An overflow error has been catched. Aborting calculation.")
                # Leave loop.
                break
    # Return the approximation of Archimedes constant.
    return ac

# Start values for the calculation of the inner polygon.
# OA = Incircle radius
# OC = Circumcircle radius
# AC = Half of edge length
# Start values 6-gon (hexagon)
# OA = 1 -> unit circle.
OA = 1
OC = OA*(2/3)*math.sqrt(3)
AC = OA/(math.sqrt(3))

# Run a simple test.
iteration = 1020
Pi = archimedes_outer_polygon(OA, OC, AC, iteration=iteration)
print(Pi)
