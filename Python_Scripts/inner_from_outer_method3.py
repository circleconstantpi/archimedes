#!/usr/bin/python3
'''Third modified Archimedes algorithm for calculating the perimeter
of the inner regular polygon based on Archimedes approach for the outer
perimeter.

Description:
Pi is calculated for the inner perimeter using Archimedes approach for
the outer perimeter of a regular polygon. Strictly spoken, this value
of Pi is a lower bound of Pi, since only the inner polygon is considered.

Neither the accuracy of the calculation with regard to the decimal
places nor the accuracy of the approximation for the circle number Pi
with a high or increasing number of iterations is taken into account.

On the test system Pi is reached with a precision of 13 places after 21
iterations.

Limitation:
So fare Pi can be calculated up to 13 correct places. The code is
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

# pylint: disable=redefined-outer-name
# pylint: disable=invalid-name

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C), 2023 Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.2"

# Initialise the number of iterations.
ITERATION = 1021

# Import the standard Python module math.
import math

# Define the function for the iterative calculation of Pi.
def inner_from_outer_method3(OB, OE, BE, iteration=5):
    '''Second modified Archimedes algorithm for calculating the perimeter
    of the inner regular polygon.'''
    # Set some strings.
    errmsg = "Aborting calculation."
    errmsg0 = "The calculation runs out of the valid range of values."
    errmsg1 = "An overflow error has occurred during the calculation."
    # Initialise the local variable.
    oldac = 0
    # Initialise the value of OA for the first iteration step.
    OA = OE
    # Run a for loop in the range from 0 to the value of iteration.
    for i in range(0, iteration+1):
        # Calculate the number of edges.
        n = 6 * 2**i
        # No iteration on first loop.
        if i == 0:
            # Calculate the approximation for Pi.
            ac = BE*n
            # Store value of pi in oldpi.
            oldac = ac
        else:
            try:
                # Calculate the length of the hypotenuse and the length of the edge.
                BF = (BE*OB)/(OB+OE)
                OF = math.sqrt(OB**2+BF**2)
                # Store the values for the next iteration.
                OB = OB+((OA-OF)*math.sqrt(1-(BF**2/OF**2)))
                OE = OA
                BE = BF*OA/OF
                # Calculate the approximation for Pi.
                ac = BE*n
                # Check whether the value for Pi is increasing.
                if ac < oldac:
                    # Print an error message.
                    print(errmsg0 + "\u0020" + errmsg)
                    # Restore value of ac.
                    ac = oldac
                    # Leave loop.
                    break
                # Store value of pi in oldpi.
                oldac = ac
            except OverflowError:
                # Print an error message.
                print(errmsg1 + "\u0020" + errmsg)
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
    # OB = Incircle radius
    # OE = Circumcircle radius
    # BE = Half of edge length
    OB = 1/2*math.sqrt(3)
    OE = 1
    BE = 1/2
    # Run a simple test.
    Pi = inner_from_outer_method3(OB, OE, BE, iteration=iteration)
    print(Pi)
    # End of function. Return None.
    return None

# Execute the script as module or as program.
if __name__ == '__main__':
    # Call the main script function.
    main(ITERATION)
