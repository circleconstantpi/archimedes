#!/usr/bin/python3
'''Applied Archimedes algorithm for calculating the perimeter of the inner
regular polygon based on Archimedes approach.

Description:
In principle Pi can be calculated up to a unknown value of correct places.
Strictly spoken, the algorithm results in a lower bound of pi, since only
the inner polygon is considered.

Limitation:
No limitations known yet.

Test system:
Python 3.8.10; pylint 2.4.4; Linux Mint 20.3 Una, Ubuntu Focal, GNU/Linux,
x86_64

Exemplary test result:
With a precision of 202 the value of Pi is calculated to 100 places
after 165 iterations.

To-Do:
Optimisation of code.

Differences to SageMath:
1. The standard Python module decimal is used.
2. Python's exponentiation operator is ** not ^.

See also:
docs.python.org/3.8/library/decimal.html
The Works of Archimedes, Measurement of a Circle
'''
# pylint: disable=useless-return
# pylint: disable=invalid-name

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023 Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.3"

# Import the standard Python module math.
from decimal import Decimal as D
from decimal import getcontext

# Initialise the constants.
ITERATION = 165
PRECISION = 202

# Set the precision of the decimal calculation.
getcontext().prec = PRECISION

# ------------------------------------------------------------------------------
# Function archimedes_inner_polygon()
# ------------------------------------------------------------------------------
def archimedes_inner_polygon(AB, AC, BC, iteration=5):
    '''Archimedes algorithm for calculating the perimeter
    of the inner regular polygon.'''
    # Run a for loop in the range from 0 to the value of iteration.
    for i in range(0, iteration+1):
        # Calculate the number of edges.
        n = 6*2**i
        # No iteration on first loop.
        if i == 0:
            # Calculate the approximation for pi.
            ac = (D(BC)*n)/2
            # Store value of pi.
            oldac = ac
        else:
            # Calculate the length of the hypotenuse and the length of the edge.
            AD = D(AB)/(D(BC**2/((AB+AC)**2)+1).sqrt())
            BD = D(AB**2 - AD**2).sqrt()
            # Store the values for the next iteration.
            BC = D(BD)
            AC = D(AD)
            # Calculate the approximation for pi.
            ac = (D(BD)*n)/2
            # The approximation of pi must be an increasing value.
            if ac <= oldac:
                # Print an error message.
                print("The value of Pi is out of range. Calculation aborted.")
                # Recover value of Pi.
                ac = oldac
                # Leave loop.
                break
            # Store value of pi.
            oldac = ac
    # Return the approximation of Archimedes constant.
    return ac

# ------------------------------------------------------------------------------
# Function correct_digits()
# ------------------------------------------------------------------------------
def correct_digits(chkpi, refpi=''):
    '''Compare a given value of Pi with a known exact reference value of Pi.
    Calculate the correct digits of a given pi number and returns a truncated
    Pi value.'''
    # Define three parts of circle number pi.
    a = "3."
    b = "14159265358979323846264338327950288419716939937510"
    c = "58209749445923078164062862089986280348253421170679"
    # Join reference pi.
    if refpi == '':
        refpi = a + b + c
    # Initialise the local variables.
    correct = ''
    idx = 0
    # Try to compare two strings.
    try:
        for char in str(chkpi):
            # Check a char against a reference char.
            if char == str(refpi)[idx]:
                # Assemble Pi.
                correct += char
                # Increment counter.
                idx += 1
            else:
                # Leave loop.
                break
    except IndexError:
        # Do nothing on error.
        pass
    # Return correct places and number.
    return (correct, idx-2)

# Main script function.
def main(precision, iteration):
    '''Main script function.'''
    # Start values 6-gon (hexagon) for the calculation of the inner polygon.
    # OB = Incircle radius
    # OE = Circumcircle radius
    # BE = Half of edge length
    AC = D(3).sqrt()
    AB = D(2)
    BC = D(1)
    # Run a simple test.
    Pi = archimedes_inner_polygon(AB, AC, BC, iteration=iteration)
    # Calculate correct number of places.
    #correct, digits = correct_digits(Pi, refpi=Pi100)
    correct, digits = correct_digits(Pi, refpi='')
    # Print summary.
    print("Precision: {}".format(precision))
    print("Iterations: {}".format(iteration))
    print("\nCalculated value of Pi rounded to 101 places:")
    print("{:.101f}".format(Pi))
    print("\nCorrect places: " + str(digits) + " ->  Pi: " + str(correct))
    # End of function. Return None.
    return None

# Execute script as module or as program.
if __name__ == '__main__':
    # Call main script function.
    main(PRECISION, ITERATION)
