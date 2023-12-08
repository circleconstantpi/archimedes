#!/usr/bin/python3
'''Zeroth modified Archimedes algorithm for calculating the perimeter
of the inner regular polygon based on Archimedes approach for the outer
perimeter.

Description:
Pi is calculated for the inner perimeter using Archimedes approach for
the outer perimeter of a regular polygon. Strictly spoken, this value
of Pi is a lower bound of Pi, since only the inner polygon is considered.

Neither the accuracy of the calculation with regard to the decimal
places nor the accuracy of the approximation for the circle number Pi
with a high or increasing number of iterations is taken into account.

On the test system Pi is reached with a precision of 14 places after 23
iterations.

Limitation:
So fare Pi can be calculated up to 14 correct places. The code is
limited to round about 1021 iterations on the test system. Then an
overflow error occurs.

To-Do:
Write a version of this script using the standard Python module
'decimal' to overcome the limitation.

Test system:
Python 3.8.10; Linux Mint 20.3 Una, Ubuntu Focal, GNU/Linux, x86_64

Differences to the usage in SageMath:
1. The standard Python module math has to be imported.
2. Python's exponentiation operator is ** not ^.
'''
# pylint: disable=invalid-name
# pylint: disable=useless-return

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023, Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.1"

# Import the standard Python module math.
import math

# Initialise the number of iterations.
ITERATION = 1021

# ----------------------------------------------------------------------
# Define the function inner_from_outer_method0()
# ----------------------------------------------------------------------
def inner_from_outer_method0(OA, OC, AC, iteration=5):
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
            ac = (1/2)*n
        else:
            # Catch an overflow error. Returns the last valid value.
            try:
                # Calculate the length of the hypotenuse and the length of the edge.
                AD = (AC*OA)/(OA+OC)
                OD = math.sqrt(OA**2+AD**2)
                # Store the values for the next iteration.
                AC = AD
                OC = OD
                # Calculate the approximation for pi.
                ac = (AD/OD/r)*n
            except OverflowError:
                # Print an error message.
                print("An overflow error has been catched. Aborting calculation.")
                # Leave loop.
                break
    # Return the approximation of Archimedes constant.
    return ac

# ++++++++++++++++++++
# Main script function
# ++++++++++++++++++++
def main(iteration):
    '''Main script function.'''
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
    Pi = inner_from_outer_method0(OA, OC, AC, iteration=iteration)
    print(Pi)
    # End of function. Return None.
    return None

# Execute the script as module or as program.
if __name__ == '__main__':
    # Call the main script function.
    main(ITERATION)
