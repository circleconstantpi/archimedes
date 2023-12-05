#!/usr/bin/python3
'''Archimedes algorithm for calculating the perimeter of the inner regular
polygon based on Archimedes approach.

Description:
Pi can be calculated up to 6 correct places. Strictly spoken, this is a lower
bound of pi, since only the inner polygon is considered.

Limitation:
The code is limited to round about 20 iterations on the test system. 6 places
for Pi is a worst result. After 15 iterations the results run out of range.
Then Pi will become 0.0, which is definitely wrong.

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
def archimedes_inner_polygon(AB, AC, BC, iteration=5, verbose=False):
    '''Archimedes algorithm for calculating the perimeter
    of the inner regular polygon.'''
    # Run a for loop in the range from 0 to the value of iteration.
    for i in range(0, iteration+1):
        # Calculate the number of edges.
        n = 6 * 2**i
        if verbose is True:
            print("Iteration: ", i)
            print("Number of edges: ", n)
        # No iteration on first loop.
        if i == 0:
            # Calculate the approximation for pi.
            ac = (BC * n)/2
        else:
            # Calculate the length of the hypotenuse and the length of the edge.
            AD = AB/math.sqrt((BC**2/(AB + AC)**2) + 1)
            BD = math.sqrt(AB**2 - AD**2)
            if verbose is True:
                print("AD: ", AD)
                print("BD: ", BD)
            # Store the values for the next iteration.
            BC = BD
            AC = AD
            # Calculate the approximation for pi.
            ac = (BD * n)/2
    # Return the approximation of Archimedes constant.
    return ac

# Print table.
def print_table():
    '''Print table.'''
    # Start values 6-gon (hexagon) for the calculation of the inner polygon.
    AC = math.sqrt(3)
    AB = 2
    BC = 1
    print("{0:<10s} | {1:<5s} | {2:<18s}".format("Iterations", "Edges", "Pi"))
    print("{0}".format(39*"-"))
    for i in range (0,5):
        n = 6 * 2**i
        Pi = archimedes_inner_polygon(AB, AC, BC, iteration=i)
        print("{0:<10d} | {1:<5d} | {2:<.15f}".format(i, n, Pi))
    # End of function. Return None.
    return None

# Calculate Pi.
def calculate_pi():
    '''Calculate Pi.'''
    # Get user input.
    iteration = int(input("\nInput number of iterations followed by [ENTER]:\u0020"))
    # Start values 6-gon (hexagon) for the calculation of the inner polygon.
    # OB = Incircle radius
    # OE = Circumcircle radius
    # BE = Half of edge length
    AC = math.sqrt(3)
    AB = 2
    BC = 1
    # Run a simple test.
    Pi = archimedes_inner_polygon(AB, AC, BC, iteration=iteration)
    print("\n{0:<.15f}".format(Pi))
    # Return None.
    return None

# Main script function.
def main():
    '''Main script function.'''
    # Call function.
    print_table()
    # Call function.
    calculate_pi()
    # End of function. Return None.
    return None

# Execute script as module or as program.
if __name__ == '__main__':
    # Call main script function.
    main()
