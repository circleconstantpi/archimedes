#!/usr/bin/python3
'''Archimedes algorithm for calculating the perimeter of the inner regular
polygon based on Archimedes approach.

Description:
Pi can be calculated up to a unknown value of correct places. Strictly
spoken, this is a lower bound of pi, since only the inner polygon is
considered.

Limitation:
No limitations known yet.

Test system:
Python 3.8.10; Linux Mint 20.3 Una, Ubuntu Focal, GNU/Linux, x86_64

Result:
With a precision of 202 the value of Pi is calculated to 100 places
after 165 iterations.

To-Do:
Optimisation of code.

Differences to SageMath:
1. The module standard Python module decimal is used.
2. Python's exponentiation operator is ** not ^.
'''
# pylint: disable=redefined-outer-name
# pylint: disable=invalid-name

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023 Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.1"

# Import the standard Python module math.
from decimal import getcontext
from decimal import Decimal as D

# Set the precision of the decimal calculation.
precision = 202
getcontext().prec = precision

# Define the function for the iterative calculation of Pi.
def archimedes_inner_polygon(AB, AC, BC, iteration=5, verbose=True):
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
            #BD = D(D(AB)**2 - D(AD)**2).sqrt()
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

# Calculate the correct digits of a given pi number.
def correct_digits(chkpi, refpi=''):
    '''Calculate the correct digits of a given pi number.'''
    # Define three parts of circle number pi.
    a = "3."
    b = "14159265358979323846264338327950288419716939937510"
    c = "58209749445923078164062862089986280348253421170679"
    if refpi == '':
        refpi = a + b + c
    correct = ''
    idx = 0
    try:
        for char in str(chkpi):
            if char == str(refpi)[idx]:
                correct += char
                idx += 1
            else:
                break
    except:
        pass
    return (correct, idx-2)

# Define Pi with 100 places.
Pi100= "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"

# Start values 6-gon (hexagon) for the calculation of the inner polygon.
# OB = Incircle radius
# OE = Circumcircle radius
# BE = Half of edge length
AC = D(3).sqrt()
AB = D(2)
BC = D(1)

# Run a simple test.
iteration = 165
Pi = archimedes_inner_polygon(AB, AC, BC, iteration=iteration)

# Calculate correct number of places.
correct, digits = correct_digits(Pi, refpi=Pi100)

# Print summary.
print("Precision:", precision)
print("Iterations:", iteration)
print("\nCalculated value of Pi rounded to 101 places:")
print("{:.101f}".format(Pi))
print("\nCorrect places:", digits, " ->  Pi:", correct)
