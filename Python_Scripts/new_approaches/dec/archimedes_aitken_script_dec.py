#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''Archimedes iterative algorithm using the Aitken's delta-squared
process.

Limitations:
If the precision is too low, errors can be triggered when AX values
should be calculated with Aitken's delta-squared process.

See also:
en.wikipedia.org/wiki/Aitken's_delta-squared_process
encyclopediaofmath.org/wiki/Aitken_Delta^2_process
'''
# pylint: disable=invalid-name
# pylint: disable=global-statement
# pylint: disable=unused-import

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023, Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.1"

# From standard Python module import some names.
from decimal import Decimal as D
from decimal import getcontext, setcontext, Context, \
                    InvalidOperation, DivisionByZero, \
                    ROUND_HALF_DOWN, ROUND_HALF_EVEN, \
                    ROUND_DOWN, ROUND_FLOOR

# Set the global constants.
ITERATION = 84
PRECISION = 152
ROUNDING = ROUND_HALF_DOWN

# Define and set the user defined context.
local_context = Context(prec=PRECISION, rounding=ROUNDING)
setcontext(local_context)

# Define a heredoc consisting of Pi with 100 places.
PI100 = '''
3.
14159265358979323846264338327950288419716939937510
58209749445923078164062862089986280348253421170679
'''

# ----------------------------------------------------------------------
# Helper function remove_ws()
# ----------------------------------------------------------------------
def remove_ws(instr):
    '''Remove whitespaces defined by a list.
    '''
    # Define the whitespaces to remove.
    mapping = [("\n", ""), ("\r", ""), ("\t", ""), (" ", "")]
    # Remove the whitespaces.
    for k, v in mapping:
        instr = instr.replace(k, v)
    # Return trimmed string.
    return instr

# ----------------------------------------------------------------------
# Function archimedes_aitken()
# ----------------------------------------------------------------------
def archimedes_aitken(iteration):
    '''Calculate Pi using the Archimedes algorithm and the Aitken's
    delta-squared process
    '''
    # Define the start values.
    a0 = D(2) * D(3).sqrt()  # half of outer perimeter
    b0 = D(3)                # half of inner perimeter
    # Define the lists for saving the lower and the upper bounds.
    lower = []
    upper = []
    # Define Aitken's lambda function.
    AX = lambda x: D(x[0]*x[2] - x[1]**2) / D(x[0] + x[2] - 2*x[1])
    # Run an iteration from 0 to iteration plus 1.
    for i in range(0, iteration+1):
        # No calculation on first run.
        if i == 0:
            a1 = D(a0)
            b1 = D(b0)
        else:
            a1 = D(2*a0*b0)/D(a0 + b0)
            b1 = D(b0*a1).sqrt()
        # Store the calculated values for the next iteration.
        b0 = D(b1)
        a0 = D(a1)
        # Add lower and upper bounds to the given lists.
        lower.append(b1)
        upper.append(a1)
        # Aitken can be used up from 3 list elements.
        if i >= 2:
            # Try to calculate Aitken's AX values.
            try:
                a2 = AX(upper)
                b2 = AX(lower)
            except InvalidOperation:
                # On error ignore Aitken's AX values.
                a2, b2 = D(a1), D(b1)
            except DivisionByZero:
                # On error ignore Aitken's AX values.
                a2, b2 = D(a1), D(b1)
            # Remove the first list element.
            lower.pop(0)
            upper.pop(0)
        else:
            # Do not calculate Aitken's AX values.
            a2, b2 = D(a1), D(b1)
    # Calculate Archimedes' constant using the arithmetic mean.
    ac = D(b2 + a2) / D(2)
    # Return Archimedes' constant.
    return ac

# ++++++++++++++++++++
# Main script function
# ++++++++++++++++++++
def main(iteration):
    '''Main script function.'''
    # Declare the global variable.
    global PI100
    # Remove whitespaces from herestring.
    PI100 = remove_ws(PI100)
    # Call function.
    ac = archimedes_aitken(iteration)
    # Print result to screen.
    print("Calc:", str(ac)[:102])
    print("Ref: ", PI100)
    # End of function. Return 1.
    return 1

# Execute script as module or as program.
if __name__ == '__main__':
    # Call main script function.
    main(ITERATION)
