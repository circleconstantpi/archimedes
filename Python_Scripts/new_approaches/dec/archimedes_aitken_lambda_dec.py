#!/usr/bin/python3
'''Archimedes iterative algorithm using the Aitken's delta-squared
process.
'''

# Import the standard Python module math.
from decimal import Decimal as D
from decimal import getcontext

# Set the constants.
PRECISION = 152
ITERATION = 84

# Set the precision of the calculation.
getcontext().prec = PRECISION

PI100 = '''
3.
14159265358979323846264338327950288419716939937510
58209749445923078164062862089986280348253421170679
'''

def remove_ws(mystr):
    mapping = [("\n", " "), ("\r", " "), (" ", "")]
    for k, v in mapping:
        mystr = mystr.replace(k, v)
    return mystr

def archimedes_aitken(iteration):
    # Define the start values.
    a0 = D(2) * D(3).sqrt()  # half of outer perimeter
    b0 = D(3)                # half of inner perimeter
    # Define to lists for saving lower and upper bound.
    lower = []
    upper = []
    # Def Aitken's lambda function.
    aitken = lambda x: D(x[0]*x[2] - x[1]**2) / D(x[0] + x[2] - 2*x[1])
    # Run an iteration from 0 to ITERATION plus 1.
    for i in range(0, ITERATION+1):
        if i == 0:
            a1 = D(a0)
            b1 = D(b0)
        else:
            a1 = D(2*a0*b0)/D(a0 + b0)
            b1 = D(b0*a1).sqrt()
        b0 = D(b1)
        a0 = D(a1)
        lower.append(b1)
        upper.append(a1)
        # Aitken can be used up from 3 list elements.
        if i >= 2:
            a2 = aitken(upper)
            b2 = aitken(lower)
            # Remove the first list element.
            lower.pop(0)
            upper.pop(0)
        else:
            a2 = D(a1)
            b2 = D(b1)
    # Calculate Archimedes' constant using the arithmetic mean.
    ac = D(b2 + a2) / D(2)
    # Return Archimedes' constant.
    return ac

# Main script function.
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
