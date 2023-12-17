#!/usr/bin/python3
'''Cubic root algorithms.'''

# Import the standard Python modules.
from decimal import Decimal as D
from decimal import getcontext, ROUND_HALF_EVEN

# Initialise the constant.
PRECISION = 64

# Set precision and rounding.
getcontext().prec = PRECISION
getcontext().rounding = ROUND_HALF_EVEN

# ----------------------------------------------------------------------
# Function cubic_root_v0()
# ----------------------------------------------------------------------
def cubic_root_v0(decnum):
    '''Applying the Newton-Raphson method to get the cubic root.'''
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec
    # Set the convergence criterion.
    eps = D(10)**(-(D(prec)-1))
    # Initialise the local variables.
    d2 = D(2)
    d3 = D(3)
    # Set and calculate the start values.
    x0 = decnum
    xn = (d2*x0 + decnum/D(x0*x0)) / d3
    # Iterate until convergence condition is reached.
    while abs(xn-x0) > eps:
        x0 = xn
        xn = (d2*x0 + decnum/D(x0*x0)) / d3
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v1()
# ----------------------------------------------------------------------
def cubic_root_v1(decnum):
    '''Applying the Newton-Raphson method to get the cubic root.'''
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec
    # Set the convergence criterion
    eps = D(10)**(-(D(prec)-1))
    # Initialise the local variables.
    d1 = D(1)
    d2 = D(2)
    d3 = D(3)
    # Calculate the start values.
    x0 = (decnum - d1)/D(d3)
    xn = (d2*x0 + decnum/D(x0*x0)) / D(d3)
    # Run the iteration until the exit condition is reached.
    while abs(xn-x0) > eps:
        x0 = xn
        xn = (d2*x0 + decnum/D(x0*x0)) / d3
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v2()
# ----------------------------------------------------------------------
def cubic_root_v2(decnum):
    '''Use binary search for determining the cubic root.'''
    # Define the helper function.
    def diff(num, mid):
        if num > (mid * mid * mid):
            return (num - (mid * mid * mid))
        else:
            return ((mid * mid * mid) - num)
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec
    # Set the convergence criterion.
    eps = D(10)**(-(D(prec)-2))
    # Set low and high for the binary search.
    low = D(0)
    high = D(decnum)
    # Run an infinite loop.
    while True:
        # Calculate the value in the middle.
        mid = low + (high - low) / 2
        # Calculate the error value.
        error = diff(decnum, mid)
        # If the value of error is less than eps leave loop.
        if error <= eps:
            break
        # Check the value of mid*mid*mid against decnum.
        if (mid * mid * mid) > decnum:
            high = mid
        else:
            low = mid
    # Return the cubic root.
    return mid

# ----------------------------------------------------------------------
# Function cubic_root_v3()
# ----------------------------------------------------------------------
def cubic_root_v3(num):
    '''Use bisection for determining the cubic root.'''
    # Define the governing function.
    def f0(x, num):
        return (x*x*x - num)
    # Define the derivative of the governing function.
    def f1(x):
        return (3 * x*x)
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec
    # Set the convergence criterion.
    eps = D(10)**(-(D(prec)-2))
    # Define the range within the result can be found
    low = D(-num)
    high = D(num)
    # Set the start value.
    x = D(0)
    # Run an infinite loop.
    while True:
        # Calculate the average value.
        x = (low + high) / 2
        # Calculate function value as well as the function derivative value.
        fvalue = f0(x, num)
        dvalue = f1(x)
        # Calculate a new range.
        if fvalue * dvalue <= 0:
            low = x
        else:
            high = x
        # Leave loop if cubic root is found.
        if fvalue < eps and fvalue >= 0:
            break
    # Return the cubic root.
    return x


# Test values.
a = D(3.46410161513775)
b = D(3.21539030917347)

# Calculate the mean value.
c = D(a + b) / D(2)

# Exponentiate c.
d = D(c)*D(c)*D(c)

# Print c.
print(c)

# Print the results of the cbic root.
print(cubic_root_v0(d))
print(cubic_root_v1(d))
print(cubic_root_v2(d))
print(cubic_root_v3(d))
