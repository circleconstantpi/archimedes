#!/usr/bin/python3
'''Third modified Archimedes algorithm for calculating the perimeter
of the inner regular polygon based on Archimedes approach for the outer
perimeter.

Description:
Pi is calculated for the inner perimeter using Archimedes approach for
the outer perimeter of a regular polygon. Strictly spoken, this value
of Pi is a lower bound of Pi, since only the inner polygon is considered.

A unit circle as well as a circle with arbitrary radius can be used for
the calculation.

On the test system Pi is reached with a precision of 13 places after 21
iterations.

Limitation:
So fare Pi can be calculated up to 13 correct places. The code is
limited to round about 1021 iterations on the test system. Then an
overflow error occurs. Decimal places are limited to 15-17 decimal
places due to the specification of floats.

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
# pylint: disable=too-many-locals
# pylint: disable=too-many-arguments

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C), 2023 Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.4"

# Import the standard Python module math.
import math

# Initialise the number of iterations.
ITERATION = 1021

# Set the ouput flags.
VERBOSE = True
SILENT = False

# ******************************************
# Define the helper function verbose_print()
# ******************************************
if VERBOSE:
    def verbose_print(*args):
        # Print each argument on a separate line.
        for arg in args:
            print(arg)
        # End of function. Return None.
        return None
else:
    # Do nothing lambda function.
    verbose_print = lambda *args: None

# *****************************************
# Define the helper function silent_print()
# *****************************************
if not SILENT:
    def silent_print(*args):
        # Print each argument on a separate line.
        for arg in args:
            print(arg)
        # End of function. Return None.
        return None
else:
    # Do nothing lambda function.
    silent_print = lambda *args: None

# ------------------------------------------------------------------------------
# Define the function inner_from_outer_method3()
# ------------------------------------------------------------------------------
def inner_from_outer_method3(OB, OE, BE, iteration=5):
    '''Third modified Archimedes algorithm for calculating the perimeter
    of the inner regular polygon.

    Arguments:
        OE (float): Circumcircle radius
        OB (float): Incircle radius
        BE (float): Half of edge length
        iteration (int): number of iterations

    Returns:
        ac (float): approximation of Pi

    OE, OB and BE are the start values for the calculation of the inner polygon.
    '''
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
        n = 6*2**i
        # Print output to the terminal.
        verbose_print("-"*24, "Loop: {}".format(i), "Edges: {}".format(n))
        # No iteration on first loop.
        if i == 0:
            # Calculate the approximation for Pi.
            ac = (BE/OA)*n
            # Store value of pi in oldpi.
            oldac = ac
            # Print output to the terminal.
            verbose_print("BE: {}".format(BE), "OA: {}".format(OA), "Pi: {}".format(ac))
        else:
            try:
                # Calculate the length of the hypotenuse and the length of the edge.
                BF = (BE*OB)/(OB+OE)
                OF = math.sqrt(OB**2+BF**2)
                # Print output to the terminal.
                verbose_print("BF: {}".format(BF), "OF: {}".format(OF))
                # Store the values for the next iteration.
                OB = OB+((OA-OF)*math.sqrt(1-(BF**2/OF**2)))
                BE = BF*OA/OF
                # Print output to the terminal.
                verbose_print("OB: {}".format(OB), "BE: {}".format(BE))
                # Calculate the approximation for Pi.
                ac = (BE/OA)*n
                # Print output to the terminal.
                verbose_print("Pi: {}".format(ac))
                # Check whether the value for Pi is increasing.
                if ac < oldac:
                    # Print output to the terminal.
                    silent_print("{0} {1}".format(errmsg0, errmsg))
                    # Restore the value of ac.
                    ac = oldac
                    # Leave loop.
                    break
                # Store value of pi in oldpi.
                oldac = ac
            except OverflowError:
                # Print output to the terminal.
                silent_print("{0} {1}".format(errmsg1, errmsg))
                # Leave loop.
                break
    # Return the approximation of the Archimedes constant.
    return ac

# ++++++++++++++++++++
# Main script function
# ++++++++++++++++++++
def main(iteration):
    '''Main script function.'''
    # Start values for the calculation of the inner polygon.
    # OE = Circumcircle radius
    # OB = Incircle radius
    # BE = Half of edge length
    # OE = 1 <- unit circle
    OE = 1
    OB = OE*(1/2)*math.sqrt(3)
    BE = OE*(1/2)
    # Run a simple test.
    Pi = inner_from_outer_method3(OB, OE, BE, iteration=iteration)
    print(Pi)
    # End of function. Return None.
    return None

# Execute the script as module or as program.
if __name__ == '__main__':
    # Call the main script function.
    main(ITERATION)
