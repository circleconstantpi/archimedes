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
Check of the code w.r.t errors. Optimisation of the code itself.
Remove typing errors.

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
__version__ = "0.1"

# Import the standard Python module math.
from decimal import Decimal as D
from decimal import getcontext

# Initialise the constants.
ITERATION = 59
PRECISION = 86

# Set the precision of the decimal calculation.
getcontext().prec = PRECISION

# Import the standard Python module math.
import math

# Define the function for the iterative calculation of Pi.
def archimedes_inner_polygon(AB, AC, BC, iteration=5, verbose=True):
    '''Archimedes algorithm for calculating the perimeter
    of inner and outer regular polygon.'''
    # Run a for loop in the range from 0 to the value of iteration.
    for i in range(0, iteration+1):
        # Calculate the number of edges.
        n = 6*2**i
        # No iteration on first loop.
        if i == 0:
            # Calculate the approximation for pi.
            ac0 = (D(BC)*n)/2
            ac1 = ((D(2)/D(3))*D(3).sqrt()*n)/2
            ac = (D(ac0)+D(ac1))/2
        else:
            # Calculate the length of the hypotenuse and the length of the edge.
            AD = D(AB)/D((D(BC)**2/D(D(AB)+D(AC))**2)+1).sqrt()
            BD = D(D(AB)**2-D(AD)**2).sqrt()
            # Store the values for the next iteration.
            BC = D(BD)
            AC = D(AD)
            # Calculate the approximation for pi.
            ac0 = (D(BD)*n)/2
            ac1 = (D(BD)/D(AD))*n
            ac = (D(ac0)+D(ac1))/2
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
Pi = archimedes_inner_polygon(AB, AC, BC, iteration=ITERATION)
print("Precision:", PRECISION)
print("Iterations:", ITERATION)
print("Ludolph van Ceulen (* 1540; † 1610) calculated 35 places. We do so, too.")
print("This calculation is like a simulation of the calculation by a human being.")
print("Calculation: {:.35f}".format(Pi))
print("Reference:   3.14159265358979323846264338327950288")
