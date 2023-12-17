#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''Cubic root algorithms.

Implementation of various procedures for study purposes,
namely Newton-Raphson, Bisection and Binary search.
'''
# pylint: disable=invalid-name
# pylint: disable=redefined-outer-name
# pylint: disable=no-else-return
# pylint: disable=chained-comparison

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023, Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.3"

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
    prec = -D(c.prec-1)
    # Set the convergence criterion.
    eps = D(10)**(prec)
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
    prec = -D(c.prec-1)
    # Set the convergence criterion
    eps = D(10)**(prec)
    # Initialise the local variables.
    d1 = D(1)
    d2 = D(2)
    d3 = D(3)
    # Calculate the start values.
    x0 = D(1)/D(3) if decnum == 1 else (decnum - d1)/D(d3)
    #x0 = (decnum - d1)/D(d3)
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
            return num - (mid * mid * mid)
        else:
            return (mid * mid * mid) - num
    # Get the used decimal precision.
    c = getcontext()
    prec = -D(c.prec-1)
    # Set the convergence criterion.
    eps = D(10)**(prec)
    # Set low and high for the binary search.
    low = D(0)
    high = D(1) if decnum < 1 else D(decnum)
    # Run an infinite loop.
    while True:
        # Calculate the value in the middle.
        mid = (low + high) / 2
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
        return x*x*x - num
    # Define the derivative of the governing function.
    def f1(x):
        return 3 * x*x
    # Get the used decimal precision.
    c = getcontext()
    prec = -D(c.prec-2)
    # Set the convergence criterion.
    eps = D(10)**(prec)
    # Define the range within the result can be found.
    low = D(0)
    high = D(1) if num < 1 else D(num)
    # Run an infinite loop.
    while True:
        # Calculate the average value.
        x = D(low + high) / 2
        # Calculate function value as well as the function derivative value.
        fvalue = f0(x, num)
        dvalue = f1(x)
        # Calculate a new range.
        if D(fvalue) * D(dvalue) <= 0:
            low = D(x)
        else:
            high = D(x)
        # Leave loop if cubic root is found.
        if (fvalue <= eps) and (fvalue >= 0):
            break
    # Return the cubic root.
    return x


# ----------------------------------------------------------------------
# Function cubic_root_v4()
# ----------------------------------------------------------------------
def cubic_root_v4(num):
    '''Use bisection for determining the cubic root.'''
    # Get the used decimal precision.
    c = getcontext()
    prec = -D(c.prec-1)
    # Set the convergence criterion.
    eps = D(10)**(prec)
    # Define the range within the answer can be found.
    low = D(0)
    high = D(1) if num < 1 else D(num)
    # Calculate value in the middle of the given range.
    mid = D(low + high) / D(2)
    # Run an infinite loop.
    while abs(mid*mid*mid-D(num)) > eps:
        # Check mid*mid*mid against num.
        if mid*mid*mid <= num:
            low = mid
        else:
            high = mid
        # Calculate value in the middle.
        mid = D(low + high) / D(2)
    # Return the cubic root.
    return mid

# ----------------------------------------------------------------------
# Function cubic_root_v5()
# ----------------------------------------------------------------------
def cubic_root_v5(x):
    '''Babylonian cubic root.'''
    # Get the used decimal precision.
    c = getcontext()
    prec = -D(c.prec-1)
    # Set the convergence criterion.
    eps = D(10)**(prec)
    # Set the constants.
    p0 = D(3)
    p1 = p0 - 1
    # Set the start values.
    y = D(x)
    w = D(1)
    # Run a loop until exit condition.
    while abs(y - w) > eps:
        y = (p1*y + w) / p0
        w = x / y**p1
    # Return the cubic root.
    return y

# Test values.
a = D(3.46410161513775)
b = D(3.21539030917347)

# Calculate the mean value.
c = D(a + b) / D(2)

# Exponentiate c.
d = D(c)*D(c)*D(c)

# Print some results of the cubic root.
print("Check for values >= 1")
print(cubic_root_v1(d))
print(cubic_root_v2(d))
print(cubic_root_v3(d))
print(cubic_root_v4(d))
print(cubic_root_v5(d))
print()
print(cubic_root_v0(D(2)))
print(cubic_root_v1(D(2)))
print(cubic_root_v2(D(2)))
print(cubic_root_v3(D(2)))
print(cubic_root_v4(D(2)))
print(cubic_root_v5(D(2)))
print("\nCheck for values <= 1")
print(cubic_root_v0(D(1)))
print(cubic_root_v1(D(1)))
print(cubic_root_v2(D(1)))
print(cubic_root_v3(D(1)))
print(cubic_root_v4(D(1)))
print(cubic_root_v5(D(1)))
print()
print(cubic_root_v0(D(0.1)))
print(cubic_root_v1(D(0.1)))
print(cubic_root_v2(D(0.1)))
print(cubic_root_v3(D(0.1)))
print(cubic_root_v4(D(0.1)))
print(cubic_root_v5(D(0.1)))
print()
print(cubic_root_v0(D(0.9)))
print(cubic_root_v1(D(0.9)))
print(cubic_root_v2(D(0.9)))
print(cubic_root_v3(D(0.9)))
print(cubic_root_v4(D(0.9)))
print(cubic_root_v5(D(0.9)))
