#!/usr/bin/python3
'''Archimedes' algorithm for calculating the perimeter of the inner
and outer regular polygon based on the Archimedean approach. The
average value of the inner and outer perimeter is an approximation
for the circle number Pi.

Description:
In principle Pi can be calculated up to a infinite number of correct
places. Standard Python is able to calculate with maximum 15-17 places.

Limitation:
No limitations are known yet.

Test system:
Python 3.8.10; Linux Mint 20.3 Una, Ubuntu Focal, GNU/Linux, x86_64

To-Do:
Check of the code with respet to errors. Optimisation of the code
itself. Remove typing errors.

See also:
docs.python.org/3.8/tutorial/floatingpoint.html
docs.python.org/3.8/library/decimal.html
'''
# pylint: disable=redefined-outer-name
# pylint: disable=invalid-name

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023, Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.3"

# Import the standard Python module math.
from decimal import Decimal as D
from decimal import getcontext

# Initialise the constants.
ITERATION = 59
PRECISION = 39

# Set the precision of the decimal calculation.
getcontext().prec = PRECISION

# Define the function for the iterative calculation of Pi.
def archimedes_outer_polygon(OA, OC, AC, iteration=4):
    '''Archimedes algorithm for calculating the perimeter of the outer
    regular polygon.'''
    # Store incircle radius for later use.
    r = D(OA)
    # Run a for loop in the range from 0 to the value of iteration plus 1.
    for i in range(0, iteration+1):
        # Calculate the number of edges.
        n = 6*2**i
        # No calculation on first loop.
        if i == 0:
            # Set the required values for the first loop.
            AD = D(AC)
            OD = D(OC)
        else:
            # Calculate the length of the hypotenuse and the length of the edge.
            AD = D(AC*OA)/D(OA+OC)
            OD = D(OA**2+AD**2).sqrt()
            # Store the values for the next iteration.
            AC = D(AD)
            OC = D(OD)
    # Calculate the approximation for pi using lower and upper bound.
    ac = (1/r+1/OD)*(AD/2)*n
    # Return the approximation of Archimedes constant.
    return ac

# Start values for the calculation of the outer polygon.
# OA = Incircle radius
# OC = Circumcircle radius
# AC = Half of edge length
# Start values 6-gon (hexagon)
# OA = 1 -> unit circle.
OA = D(1)
OC = D(OA)*(D(2)/D(3))*D(3).sqrt()
AC = D(OA)/(D(3).sqrt())

# Run a simple test.
Pi = archimedes_outer_polygon(OA, OC, AC, iteration=ITERATION)
print("Precision:", PRECISION)
print("Iterations:", ITERATION)
print("Ludolph van Ceulen (* 1540; † 1610) calculated 35 places. We do so, too.")
print("This calculation is like a simulation of the calculation by a human being.")
print("Calculation: {:.35f}".format(Pi))
print("Reference:   3.14159265358979323846264338327950288")
