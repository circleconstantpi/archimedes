#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''Collection of cubic root algorithms.

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
__version__ = "0.4"

# Import the standard Python modules.
from timeit import default_timer as timer
from datetime import timedelta

# Import the standard Python modules.
from decimal import Decimal as D
from decimal import getcontext, localcontext, ROUND_HALF_EVEN, ROUND_HALF_DOWN

# Initialise the constant.
PRECISION = 64

# Set precision and rounding.
getcontext().prec = PRECISION
getcontext().rounding = ROUND_HALF_DOWN

# ----------------------------------------------------------------------
# Function cubic_root_v0()
# ----------------------------------------------------------------------
def cubic_root_v0(a):
    '''Applying the Newton-Raphson method to get the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = a
    xn = (2*x0 + a/D(x0*x0)) / 3
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context behaviour.
        ctx.prec += 2
        ctx.rounding=ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            xn = (2*x0 + a/D(x0*x0)) / 3
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v1()
# ----------------------------------------------------------------------
def cubic_root_v1(a):
    '''Applying the Halleys method to get the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = a
    xn = x0 * (x0*x0*x0 + 2*a) / (2 * x0*x0*x0 + a)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context behaviour.
        ctx.prec += 2
        ctx.rounding=ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            xn = x0 * (x0*x0*x0 + 2*a) / (2 * x0*x0*x0 + a)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v2()
# ----------------------------------------------------------------------
def cubic_root_v2(decnum):
    '''Applying the Newton-Raphson method to get the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(decnum).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion
    eps = D(10)**(-prec)
    # Initialise the local variables.
    d1 = D(1)
    d2 = D(2)
    d3 = D(3)
    # Calculate the start values.
    x0 = D(1)/D(3) if decnum == 1 else (decnum - d1)/D(d3)
    xn = (d2*x0 + decnum/D(x0*x0)) / D(d3)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context behaviour.
        ctx.prec += 2
        ctx.rounding=ROUND_HALF_DOWN
        # Run the iteration until the exit condition is reached.
        while abs(xn-x0) > eps:
            x0 = xn
            xn = (d2*x0 + decnum/D(x0*x0)) / d3
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v3()
# ----------------------------------------------------------------------
def cubic_root_v3(num):
    '''Use binary search for determining the cubic root.'''
    # Define the function for calculating the difference value.
    def diff(num, mid):
        if num > (mid * mid * mid):
            return num - (mid * mid * mid)
        else:
            return (mid * mid * mid) - num
    # Calculate the number of leading digits.
    cln = len(str(num).split(".")[0])
    # Get the decimal precision in use.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set low and high for the binary search.
    low = D(0)
    high = D(2) if num <= 1 else D(2*num)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context behaviour.
        ctx.prec += 2
        ctx.rounding=ROUND_HALF_DOWN
        # Run an infinite loop.
        while True:
            # Calculate the value in the middle.
            mid = (low + high) / 2
            # Calculate the error value.
            error = diff(num, mid)
            # If the value of error is less than eps leave loop.
            if error <= eps:
                break
            # Check the value of mid*mid*mid against num.
            if (mid * mid * mid) >= num:
                high = mid
            else:
                low = mid
    # Restore precision.
    mid = +mid
    # Return the cubic root.
    return mid

# ----------------------------------------------------------------------
# Function cubic_root_v4()
# ----------------------------------------------------------------------
def cubic_root_v4(num):
    '''Use bisection for determining the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(num).split(".")[0])
    # Define the governing function.
    def f0(x, num):
        return x*x*x - num
    # Define the derivative of the governing function.
    def f1(x):
        return 3 * x*x
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set the range within the answer can be found.
    low = D(0)
    high = D(1) if num < 1 else D(num)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding=ROUND_HALF_EVEN
        # Run an infinite loop.
        while True:
            # Calculate the average value.
            x = D(low + high) / 2
            # Calculate function value and function derivative value.
            fvalue = f0(x, num)
            dvalue = f1(x)
            # Calculate a new range.
            if fvalue * dvalue < 0:
                low = D(x)
            else:
                high = D(x)
            # Leave loop if cubic root is found.
            if (fvalue <= eps) and (fvalue >= 0):
                break
    # Restore precision.
    x = +x
    # Return the cubic root.
    return x

# ----------------------------------------------------------------------
# Function cubic_root_v5()
# ----------------------------------------------------------------------
def cubic_root_v5(num):
    '''Use bisection for determining the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(num).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding=ROUND_HALF_DOWN
        # Define the range within the answer can be found.
        low = D(0)
        high = D(2) if num < 1 else D(2*num)
        # Calculate value in the middle of the given range.
        mid = D(low + high) / 2
        # Run an infinite loop.
        while abs(mid*mid*mid-num) >= eps:
            # Check mid*mid*mid against num.
            if mid*mid*mid <= num:
                low = mid
            else:
                high = mid
            # Calculate value in the middle.
            mid = D(low + high) / 2
    # restore the precision.
    mid = +mid
    # Return the cubic root.
    return mid

# ----------------------------------------------------------------------
# Function cubic_root_v6()
# ----------------------------------------------------------------------
def cubic_root_v6(x):
    '''Babylonian cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(x).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set the constants.
    p0 = D(3)
    p1 = D(p0 - 1)
    # Set the start values.
    y = D(x)
    w = D(1)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding=ROUND_HALF_DOWN
        # Run a loop until convergence is reached.
        while abs(w - y) >= eps:
            y = (p1*y + w) / p0
            w = x / y**p1
    # Restore the precision.
    y = +y
    # Return the cubic root.
    return y

# ----------------------------------------------------------------------
# Function cubic_root_v7()
# ----------------------------------------------------------------------
def cubic_root_v7(a):
    '''Herons cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set the start values.
    x0 = D(a)
    #xn = D(x0) - D(x0**3-a)/D(3*x0**2)
    xn = D(x0) - D(x0*x0*x0-a)/D(3*x0*x0)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding=ROUND_HALF_DOWN
        # Run a loop until convergence is reached.
        while abs(x0 - xn) >= eps:
            x0 = xn
            #xn = D(x0) - D(x0**3-a)/D(3*x0**2)
            xn = D(x0) - D(x0*x0*x0-a)/D(3*x0*x0)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v8()
# ----------------------------------------------------------------------
def cubic_root_v8(a):
    '''Applying the Housholder method to get the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = a
    xn = x0 - D(x0**3-a)/D(3*x0**2) - D((x0**3-a)**2*6*x0)/D(2*(3*x0**2)**3)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context behaviour.
        ctx.prec += 2
        ctx.rounding=ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            xn = x0 - D(x0**3-a)/D(3*x0**2) - D((x0**3-a)**2*6*x0)/D(2*(3*x0**2)**3)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# Test values.
a = D(3.46410161513775)
b = D(3.21539030917347)

# Calculate the mean value.
c = D(a + b) / D(2)

# Exponentiate c.
d = D(c)*D(c)*D(c)

# Set up a test array.
test_array = [d, 0.1, 0.3, 0.5, 0.8, 0.999999, 1, 2, 562345.983526627, 987299826534888736635.981425, '80996655413131.092562882']
#test_array = [d]

# Perform some tests.
for i in test_array:
    print("Test value:", i)
    for j in [0, 1, 2, 3, 4, 5, 6, 7]:
        msgstr = "cubic_root_v" + str(j) + "():"
        print (msgstr)
        start = timer()
        prgstr = "cubic_root_v" + str(j) + "(D(" + str(i) + "))"
        cr = eval(prgstr)
        end = timer()
        print(cr)
        print("Elapsed time:", timedelta(seconds=end-start))
    print()




