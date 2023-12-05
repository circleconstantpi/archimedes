#!/usr/bin/python3
'''Second modified Archimedes algorithm for calculating the perimeter of the
inner regular polygon based on Archimedes approach.

Description:
Pi can be calculated up to 13 correct places. Strictly spoken, this is a lower
bound of pi, since only the inner polygon is considered.

Limitation:
The code is limited to round about 1010 iterations on the test system.

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
def archimedes_inner_polygon_method2(OB, OE, BE, iteration=5):
    '''Second modified Archimedes algorithm for calculating the perimeter
    of the inner regular polygon.'''
    # Initialise the value of OA for the first iteration step.
    OA = OE
    # Run a for loop in the range from 0 to the value of iteration.
    for i in range(0, iteration):
        # Calculate the number of edges.
        n = 6 * 2**i
        # No iteration on first loop.
        if i == 0:
            # Calculate the approximation for pi.
            ac = BE * n
        else:
            # Calculate the length of the hypotenuse and the length of the edge.
            BF = (BE*OB)/(OB + OE)
            OF = math.sqrt(OB**2 + BF**2)
            # Store the values for the next iteration.
            OB = OB + ((OA-OF)*math.sqrt(1 - (BF**2/OF**2)))
            OE = OA
            BE = BF + BF*((OA - OF)/OF)
            # Calculate the approximation for pi.
            ac = BE * n
    # Return the approximation of Archimedes constant.
    return ac

# Start values for the calculation of the inner polygon.
# OB = Incircle radius
# OE = Circumcircle radius
# BE = Half of edge length
OB = 1/2*math.sqrt(3)
OE = 1
BE = 1/2

# Run a simple test.
iteration = 1010
Pi = archimedes_inner_polygon_method2(OB, OE, BE, iteration=iteration)
print(Pi)
