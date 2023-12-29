#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''Collection of algorithms for the calculation of the cubic root.

Description:
Implementation of various procedures for the high precision calculation
of cubic roots for study purposes.

Objective:
The objective is studying which one of the algorithms is fast and
reliable. For a long term calculation where per iteration the cubic
root is needed I need a fast one. The algorithm must still be robust.

Following famous algorithms are implemented:

  Regula Falsi Method (False Position Method,  Method of False Position)
  Bisection Method
  Binary Search Method
  Newton's Method
  Newton–Raphson Method
  Halley's Method
  Householder's Method
  Secant Method
  Babylonian Method
  Heron's Method
  Steffensen Method
  Steffensen Method using Aitken squared delta method
  Brute Force Method

To-Do:

Optimisation of the code. Finding and correcting typos. Improving of
the documentation. Checking for errors. Calculate more examples.
Optimisation of code. Find and coreect typos. Improvement of the
documentation.

Take a look on following methods/methodologys:

  Bracket Approaches (e.g. Bisection Method)
  Iterative Approaches (e.g. Newton's Method)
  Fixed Point Iteration Method (e.g. Newton’s Method)
  Hybrid Methods (e.g. Brent’s Method)
  Wegstein's Method
  Illinois Method
  Pegasus Method
  Anderson-Björck Method
  Brent's Method
  Ridders' Method
  Quasi-Newton Methods
  Broyden’s Method
  ITP Method
  Bakhshali Method
  Goldschmidt’s algorithm
  Digit-by-digit calculation

Limitations:
Steffensens methods works not as expected for bigger numbers.

Usage: (from shell prompt)
    python3 ./cubic_root.py

See also:
https://www.frassek.org/numerik/nullstellen-reeller-funktionen/
www.math-cs.gordon.edu/courses/mat342/python/findroot.py
cocalc.com/share/public_paths/embed/bf8c2613e8e9f04f4b052a83fb30a14e1ccea697
'''
# pylint: disable=no-else-return
# pylint: disable=redefined-outer-name
# pylint: disable=invalid-name
# pylint: disable=chained-comparison
# pylint: disable=invalid-unary-operand-type
# pylint: disable=line-too-long
# pylint: disable=eval-used
# pylint: disable=missing-function-docstring

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023, Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.5"

# Import the standard Python modules.
from timeit import default_timer as timer
from datetime import timedelta

import decimal

# Import the standard Python modules.
from decimal import Decimal as D
from decimal import getcontext, localcontext, \
                    DivisionByZero, InvalidOperation, \
                    ROUND_HALF_EVEN, ROUND_HALF_DOWN

# Initialise the constant for the calculation precision.
PRECISION = 128

# Set precision and rounding.
getcontext().prec = PRECISION
getcontext().rounding = ROUND_HALF_DOWN

# ----------------------------------------------------------------------
# Function approximation()
# ----------------------------------------------------------------------
def approximation(a0, b0, num, prec):
    '''Calculate a approximation for the cubic root for a given radical.

    Brute force method for finding the cubic root to the number of
    requested places.

    arguments:
        num (decimal) : radical of the searched root
        a0, b0 (int)  : estimation lower and upper bound
        prec (int)    : precision

    Returns:
        a0, b0 (decimal) : lower and upper bound
    '''
    # Get the lenght of the digits before the decimal point.
    intlen = lambda c: len(str(c).split(".")[0])
    # Calculate and set the number of places.
    a0_len, b0_len = intlen(a0), intlen(b0)
    digits = prec - int(((a0_len + b0_len)/2))
    # Convert lower and upper bound to decimal.
    a0, b0 = D(a0), D(b0)
    # Initialise the initial digit.
    place = D(1)
    # Create place by place.
    for _ in range(0, digits+1):
        # Calculate the next decimal place.
        place /= 10
        # Add and subtract a new digit to a0 and b0.
        while a0*a0*a0 < num:
            a0 += place
        while b0*b0*b0 > num:
            b0 -= place
        # Decrement and increment a0 and b0 by the value of place.
        a0 -= place
        b0 += place
    # Return the lower and the upper bound.
    return a0, b0

# ----------------------------------------------------------------------
# Function search_interval()
# ----------------------------------------------------------------------
def search_interval(c):
    '''Calculate the search interval for the cubic root.

    arguments:
        c (decimal) : radical of the searched root

    Returns:
        a0, b0 (decimal) : lower and upper bound
    '''
    # Get the value of the exponent.
    intlen = len(str(c).split(".")[0])
    exp = round(((intlen-2)/3) + 1)
    # Set initial values min and max for the interval.
    minval = 0
    maxval = 10**exp
    # Initialise the return values.
    a1, b1 = int(minval), int(maxval)
    # Define a function.
    def determine_limits(c, a0, b0, step):
        '''Limit determination.'''
        # Initialise the return values.
        a1, b1 = a0, b0
        # Check if i to the power of three is larger or smaller as c.
        for i in range(a0, b0, step):
            if i*i*i < c:
                a1 = i
            elif i*i*i > c:
                b1 = i
                break
        # Return the limits of the interval.
        return a1, b1
    # Define a function.
    def reverse_search(a0, b0, maxval):
        '''Reverse search.'''
        # Initialise the return values.
        a1, b1 = a0, b0
        # Run a reverse loop until we reach the potency value 10.
        while maxval > 1:
            # Assign max to the division by 10.
            maxval /= 10
            # The calculated max value is the new step width.
            step = int(maxval)
            # Calculate the values of the interval limits.
            a1, b1 = determine_limits(c, a0, b0, step)
            a0, b0 = a1, b1
        # Return the limits of the interval.
        return a1, b1
    # Return the limits based on the value of c.
    if c > 1:
        a1, b1 = reverse_search(a1, b1, maxval)
    elif c == 1:
        a1, b1 = 0.5, 1.5
    elif c < 1:
        a1, b1 = 0, 1
    return a1, b1

# ----------------------------------------------------------------------
# Function cubic_root_v0()
# ----------------------------------------------------------------------
def cubic_root_v0(a):
    '''Applying the Newton-Raphson method to calculate the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Initialise the start values.
    x0 = D(a)
    xn = D(0)
    diff = D(1)
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until the convergence condition is fulfilled.
        while diff >= eps:
            xn = x0 - (x0*x0*x0 - a) / (3 * x0*x0)
            diff = abs(xn-x0)
            x0 = xn
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v1()
# ----------------------------------------------------------------------
def cubic_root_v1(a):
    '''Applying the Newton-Raphson method to calculate the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Initialise the start values.
    x0 = D(a)
    xn = D(0)
    diff = D(1)
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until the convergence condition is fulfilled.
        while diff >= eps:
            xn = (2 * x0*x0*x0 + a) / (3 * x0*x0)
            diff = abs(xn-x0)
            x0 = xn
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v2()
# ----------------------------------------------------------------------
def cubic_root_v2(a):
    '''Applying the Newton-Raphson method to get the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = D(a)
    xn = x0 - (x0*x0*x0 - a) / (3 * x0 * x0)
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            xn = x0 - (x0*x0*x0 - a) / (3 * x0 * x0)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v3()
# ----------------------------------------------------------------------
def cubic_root_v3(a):
    '''Applying the Newton-Raphson method to get the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = D(a)
    xn = (2*x0 + D(a)/D(x0*x0)) / 3
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            xn = (2*x0 + D(a)/D(x0*x0)) / 3
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v4()
# ----------------------------------------------------------------------
def cubic_root_v4(a):
    '''Newton-Raphson method for the calculation of a cubic root.

    Arguments:
        dec (decimal): radical of the cubic root

    Returns:
        xn (decimal): cubic root of the given radical
    '''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion
    eps = D(10)**(-prec)
    # Initialise the local variables.
    d2 = D(2)
    d3 = D(3)
    # Calculate the start values.
    x0 = D(a)
    xn = D(d2*x0 + D(a)/D(x0*x0)) / d3
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context behaviour.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Run the iteration until the exit condition is reached.
        while abs(xn-x0) > eps:
            x0 = xn
            xn = D(d2*x0 + D(a)/D(x0*x0)) / d3
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v5()
# ----------------------------------------------------------------------
def cubic_root_v5(a):
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
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            xn = x0 * (x0*x0*x0 + 2*a) / (2 * x0*x0*x0 + a)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v6()
# ----------------------------------------------------------------------
def cubic_root_v6(a):
    '''Applying the Halleys method to get the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Define the decimal number.
    d2 = D(2)
    # Set and calculate the start values.
    x0 = a
    xn = x0 * (x0*x0*x0 + 2*a) / (2 * x0*x0*x0 + a)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context behaviour.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            xn = x0 * (x0*x0*x0 + d2*a) / (d2 * x0*x0*x0 + a)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v7()
# ----------------------------------------------------------------------
def cubic_root_v7(a):
    '''Applying Halleys method for calculate the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Define the decimal number.
    d2 = D(2)
    # Set and calculate the start values.
    x0 = a
    x3 = x0*x0*x0
    xn = x0 * (x3 + d2*a) / (d2*x3 + a)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context behaviour.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            x3 = x0*x0*x0
            xn = x0 * (x3 + d2*a) / (d2*x3 + a)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v8()
# ----------------------------------------------------------------------
def cubic_root_v8(num):
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
        ctx.rounding = ROUND_HALF_DOWN
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
# Function cubic_root_v9()
# ----------------------------------------------------------------------
def cubic_root_v9(num):
    '''Use Bisection for determining the cubic root.'''
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
        ctx.rounding = ROUND_HALF_EVEN
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
# Function cubic_root_v10()
# ----------------------------------------------------------------------
def cubic_root_v10(num):
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
        ctx.rounding = ROUND_HALF_DOWN
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
# Function cubic_root_v11()
# ----------------------------------------------------------------------
def cubic_root_v11(x):
    '''Babylonian Method for calculating the cubic root.'''
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
        ctx.rounding = ROUND_HALF_DOWN
        # Run a loop until convergence is reached.
        while abs(w - y) >= eps:
            y = (p1*y + w) / p0
            w = x / y**p1
    # Restore the precision.
    y = +y
    # Return the cubic root.
    return y

# ----------------------------------------------------------------------
# Function cubic_root_v12()
# ----------------------------------------------------------------------
def cubic_root_v12(a):
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
    xn = D(x0) - D(x0*x0*x0-a)/D(3*x0*x0)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Run a loop until convergence is reached.
        while abs(x0 - xn) >= eps:
            x0 = xn
            xn = D(x0) - D(x0*x0*x0-a)/D(3*x0*x0)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v13()
# ----------------------------------------------------------------------
def cubic_root_v13(a):
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
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            xn = x0 - D(x0**3-a)/D(3*x0**2) - D((x0**3-a)**2*6*x0)/D(2*(3*x0**2)**3)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v14()
# ----------------------------------------------------------------------
def cubic_root_v14(num):
    '''Applying Regula falsi to get the cubic root.'''
    # Definition of governing equation.
    def f(x, a):
        return x*x*x - a
    # Calculate the number of leading digits.
    cln = len(str(num).split(".")[0])
    # Get the used decimal precision.
    gc = getcontext()
    prec = gc.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Initialise the variable sign.
    sign = 0
    # a, b are endpoints of an interval where we are searching.
    a = -2*num
    b = 2*num
    # Compute the function values at the endpoints of the interval.
    fa = f(a, num)
    fb = f(b, num)
    # Perform a calculation in a new block using a local context.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Run an infite loop.
        while True:
            # Calculate the estimation.
            c = (fa * b - fb * a) / (fa - fb)
            # Leave loop on convergence.
            if abs(b - a) < eps * abs(b + a):
                break
            # Calculate the function value at the midpoint of the interval
            fc = f(c, num)
            # Leave loop, when fc * fa, fb is very small or is nearly zero.
            if fc*fb > 0:
                # fc and fb have the same sign.
                # Exchange value and function value.
                b = c
                fb = fc
                # Check the sign and half the function value.
                if sign == -1:
                    fa /= 2
                sign = -1
            elif fa*fc > 0:
                # fc and fa have the same sign.
                # Exchange value and function value.
                a = c
                fa = fc
                # Check the sign and half the function value.
                if sign == +1:
                    fb /= 2
                sign = +1
            else:
                # fc * fa, fb is very small or is nearly zero.
                break
    # Restore precision.
    c = +c
    # Return the cubic root.
    return c

# ----------------------------------------------------------------------
# Function cubic_root_v15()
# ----------------------------------------------------------------------
def cubic_root_v15(num):
    '''Secant method for calculating the cube root.'''
    def f(x, num):
        return x*x*x - num
    def diff2(x):
        return 6 * x
    lb, ub = search_interval(num)
    lb, ub = approximation(lb, ub, num, 2)
    # Calculate the number of leading digits.
    cln = len(str(num).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set the borders.
    a = D(lb)
    b = D(ub)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context behaviour.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_EVEN
        # Calculate the cubic root.
        if f(b, num)*diff2(b) >= 0:
            x0 = a
            while abs(f(x0, num)) >= eps:
                x0 = x0 - (f(x0, num) / (f(b, num) - f(x0, num))) * (b - x0)
        elif f(a, num)*diff2(a) >= 0:
            x0 = b
            while abs(f(x0, num)) >= eps:
                x0 = x0 - (f(x0, num) / (f(x0, num) - f(a, num))) * (x0 - a)
    # Restore precision.
    x0 = +x0
    # Return cubic root.
    return x0

# ----------------------------------------------------------------------
# Function cubic_root_v16()
# ----------------------------------------------------------------------
def cubic_root_v16(num):
    '''Brute force method.'''
    # Get the used global decimal precision.
    c = getcontext()
    prec = c.prec
    # Calculate the search interval.
    a0, b0 = search_interval(num)
    # Run a calculation block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Calculate the best approximation for lower and upper limit.
        a0, b0 = approximation(a0, b0, num, prec)
        # Calculate the mean value of lower and upper limit.
        cbrt = D(a0 + b0) / 2
    # Restore the precision.
    cbrt = +cbrt
    # Return the cubic root.
    return cbrt

# ----------------------------------------------------------------------
# Function cubic_root_v17()
# ----------------------------------------------------------------------
def cubic_root_v17(dec):
    '''Steffensen Method for finding cubic roots.'''
    # Define the governing function f(x).
    def f(x, dec):
        fx = x*x*x - dec
        return fx
    # Define the slope function g(x)
    def g(x, dec):
        h = f(x, dec)
        gx = (f(x + h, dec) - f(x, dec)) / f(x, dec)
        return gx
    # Calculate the number of leading digits.
    cln = len(str(dec).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    dec = D(dec)
    x0 = D(1)/D(3)*dec
    x1 = x0 - f(x0, dec) / g(x0, dec)
    with localcontext() as ctx:
        # Change the local context behaviour.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        while abs((x1-x0)/x0) >= eps:
            try:
                x0 = x1
                x1 = x0 - f(x0, dec) / g(x0, dec)
            except InvalidOperation:
                pass
    # Restore the precision.
    x1 = +x1
    # Rteurn the cubic root.
    return x1

# ----------------------------------------------------------------------
# Function cubic_root_v18()
# ----------------------------------------------------------------------
def cubic_root_v18(a):
    '''Steffensen-Aitken Method for finding cubic roots.

    See also:
    en.wikipedia.org/wiki/Steffensen's_method
    www.math-cs.gordon.edu/courses/mat342/python/findroot.py
    '''
    # Define the governing equation.
    def f(x, a):
        return x*x*x - a
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Calculate the start value.
    p0 = D(1)/D(3)*D(a)
    # Run calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Run an infinite loop.
        while True:
            # Try to calculate the cubic root.
            try:
                p1 = p0 + f(p0, a)  # non-standard: adding p0 to function value
                p2 = p1 + f(p1, a)  # non-standard: adding p1 to function value
                p = p2 - (p2 - p1)*(p2 - p1) / (p2 - 2*p1 + p0)
                # Convergence criteria.
                if abs(p-p0) <= eps:
                    break
                # Store value for next loop.
                p0 = p
            except DivisionByZero:
                break
    # Restore the precision.
    p = +p
    # Rteurn the cubic root.
    return p

# ----------------------------------------------------------------------
# Function cubic_root_v19()
# ----------------------------------------------------------------------
def cubic_root_v19(a):
    '''Applying the Kou-Li-Wang method to get the cubic root.'''
    # Define the governing function.
    def f0(x, num):
        return x*x*x - num
    # Define the 1. derivative of the governing function.
    def f1(x):
        return 3 * x*x
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = a
    xn = (2*x0 + D(a)/D(x0*x0)) / 3
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            xs = x0 - (f0(x0, a) / f1(x0))
            u1 = x0 - (2*f0(x0, a)) / (f1(x0) + f1(xs))
            xn = u1 - (f0(u1, a) / f1(u1))
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v20()
# ----------------------------------------------------------------------
def cubic_root_v20(a):
    '''Applying the Weerakoon-Fernando method to get the cubic root.'''
    # Define the governing function.
    def f0(x, num):
        return x*x*x - num
    # Define the 1. derivative of the governing function.
    def f1(x):
        return 3 * x*x
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = a                           # check this for this method
    xn = (2*x0 + D(a)/D(x0*x0)) / 3  # check this for this method
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            xs = x0 - f0(x0, a) / f1(x0)
            xn = x0 - 2*f0(x0, a) / (f1(x0) + f1(xs))
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v21()
# ----------------------------------------------------------------------
def cubic_root_v21(a):
    '''Applying the Kasturiarachi method to get the cubic root.'''
    # Define the governing function.
    def f0(x, num):
        return x*x*x - num
    # Define the 1. derivative of the governing function.
    def f1(x):
        return 3 * x*x
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = a                           # check this for this method
    xn = (2*x0 + D(a)/D(x0*x0)) / 3  # check this for this method
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            # eps condition musst be modified to skip the try/ecept block.
            try:
                x0 = xn
                xs = x0 - f0(x0, a) / f1(x0)
                xn = x0 - ((f0(x0, a) * f0(x0, a)) / (f1(x0) * (f0(x0, a) - f0(xs, a))))
            except InvalidOperation:
                break
            except DivisionByZero:
                break
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v22()
# ----------------------------------------------------------------------
def cubic_root_v22(a):
    '''Applying the Schröder method to get the cubic root.'''
    # Define the governing function.
    f0 = lambda x, num: x*x*x - num
    # Define the 1. derivative of the governing function.
    f1 = lambda x: 3 * x*x
    # Define the 2. derivative of the governing function.
    f2 = lambda x: 6 * x
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = a                           # check this for this method
    xn = (2*x0 + D(a)/D(x0*x0)) / 3  # check this for this method
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            nr = f0(x0, a) * f1(x0)
            dr = f1(x0)*f1(x0) - f0(x0, a) * f2(x0)
            xn = x0 - (nr / dr)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v23()
# ----------------------------------------------------------------------
def cubic_root_v23(a):
    '''Applying the Chebyshev method to get the cubic root.'''
    # Define the governing function.
    def f0(x, num):
        return x*x*x - num
    # Define the 1. derivative of the governing function.
    def f1(x):
        return 3 * x*x
    # Define the 2. derivative of the governing function.
    def f2(x):
        return 6 * x
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = a
    xn = (2*x0 + D(a)/D(x0*x0)) / 3
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            xn = x0 - ((1 + (D(1)/D(2)) * ((f2(x0) * f0(x0, a)) / (f1(x0))**2)) * f0(x0, a)/f1(x0))
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v24()
# ----------------------------------------------------------------------
def cubic_root_v24(a):
    '''Applying the Chebyshev method to get the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = D(a)
    xn = (2*x0 + D(a)/D(x0*x0)) / 3
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            #xn = (5*x0**6 + 5*a*x0**3 - a**2)/(9*x0**5)
            xn = (5*x0*x0*x0*x0*x0*x0 + 5*x0*x0*x0*a - a**2)/(9*x0*x0*x0*x0*x0)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v25()
# ----------------------------------------------------------------------
def cubic_root_v25(a):
    '''Applying the Ostrowski method to get the cubic root.'''
    # Define the governing function.
    f0 = lambda x, num: x*x*x - num
    # Define the 1. derivative of the governing function.
    f1 = lambda x: 3 * x*x
    # Define the 2. derivative of the governing function.
    f2 = lambda x: 6 * x
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = a                           # check this for this method
    xn = (2*x0 + D(a)/D(x0*x0)) / 3  # check this for this method
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            try:
                x0 = xn
                y1 = x0 - (f0(x0, a)/f1(x0))
                nr = f0(y1, a) - f0(x0, a)
                dr = 2*f0(y1, a) - f0(x0, a)
                xn = x0 - (f0(x0, a)/f1(x0)) * (nr / dr)
            except InvalidOperation:
                break
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v26()
# ----------------------------------------------------------------------
def cubic_root_v26(a):
    '''Applying the Weerakoon-Fernando method to get the cubic root.'''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = D(a)
    xn = (2*x0 + D(a)/D(x0*x0)) / 3
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            xn = (7*x0**7 + 10*a*x0**4 + a**2*x0)/(13*x0**6 + 4*a*x0**3 + a**2)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v27()
# ----------------------------------------------------------------------
def cubic_root_v27(a):
    '''Applying the Taylor series expansion to get the cubic root.'''
    # Define factorial function.
    def factorial(n):
        if n < 2:
            return 1
        else:
            return n * factorial(n-1)
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = D(0)
    xn = D(0)
    n = 0
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while True:
            if a == 1:
                xn = 1
                break
            else:
                x0 = D((D(1)/D(3) * D(a).ln())**n) / D(factorial(n))
                xn += x0
                n += 1
            if abs(x0) <= eps:
                break
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v28()
# ----------------------------------------------------------------------
def cubic_root_v28(a):
    '''Applying an uncommon approach to get the cubic root.

    See also:
    math.stackexchange.com/questions/1400263/how-to-make-this-cubic-root-c-algorithm-faster
    '''
    # Calculate the number of the leading digits of the radical.
    cln = len(str(a).split(".")[0])
    # Get the used global decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the value for the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = D(0)
    x1 = D(a) - D(0.1)
    xn = D(a) + D(0.1)
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until the convergence condition is reached.
        while abs(xn-x0) >= eps:
            x0 = xn
            x1 = xn * (a + x1) / (x1 + xn)
            xn = x1 / xn
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v29()
# ----------------------------------------------------------------------
def cubic_root_v29(a):
    '''Cubic root following Heron's book Metrica determined form the
    square root.

    See also:
    stackoverflow.com/questions/30381552/cube-root-using-herons-method
    '''
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set the start values.
    x0 = D(a)
    xn = (D(x0) + D(a)/D(x0)) / D(2)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Run a loop until convergence is reached.
        while abs(x0 - xn) >= eps:
            x0 = xn
            xn = (3*D(x0) + D(a)/D(x0*x0)) / D(4)
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v30()
# ----------------------------------------------------------------------
def cubic_root_v30(a):
    '''Cubic root determined using Jain's method.'''
    # Define the lambda functions.
    f0 = lambda x, a: x*x*x - a
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set the start values.
    x0 = D(a)
    xn = D(a)
    diff = D(0.1)
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Run a loop until convergence is reached.
        while diff >= eps:
            try:
                h = x0 + f0(x0, a)
                fh = f0(h, a) - f0(x0, a)
                f2 = f0(x0, a)*f0(x0, a)
                f3 = f0(x0, a)*f0(x0, a)*f0(x0, a)
                xs = x0 - (D(f2) / D(fh))
                xn = x0 - (D(f3) / D(fh*(f0(x0, a) - f0(xs, a))))
                diff = abs(xn-x0)
                x0 = xn
            except InvalidOperation:
                if x0 == a and x0 != 1:
                    xn = -1
                break
            except DivisionByZero:
                if x0 == a and x0 != 1:
                    xn = -1
                break
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ----------------------------------------------------------------------
# Function cubic_root_v31()
# ----------------------------------------------------------------------
def cubic_root_v31(a):
    '''Applying the Ptak method to get the cubic root.'''
    # Define the governing function.
    def f0(x, num):
        return x*x*x - num
    # Define the 1. derivative of the governing function.
    def f1(x):
        return 3 * x*x
    # Calculate the number of leading digits.
    cln = len(str(a).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    eps = D(10)**(-prec)
    # Set and calculate the start values.
    x0 = a
    yn = x0 - f0(x0, a) / f1(x0)
    xn = x0 - (f0(x0, a) + f0(yn, a)) / f1(x0)
    # Perform the calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Iterate until convergence condition is reached.
        while abs(xn-x0) >= eps:
            # eps condition musst be modified to skip the try/ecept block.
            try:
                x0 = xn
                yn = x0 - f0(x0, a) / f1(x0)
                xn = x0 - (f0(x0, a) + f0(yn, a)) / f1(x0)
            except InvalidOperation:
                break
    # Restore the precision.
    xn = +xn
    # Return the cubic root.
    return xn

# ++++++++++++++++++++
# Main script function
# ++++++++++++++++++++
def main():
    '''Main script function.'''
    # Test values.
    a = D(3.46410161513775)
    b = D(3.21539030917347)
    # Calculate the mean value.
    c = D(a + b) / D(2)
    # Exponentiate c.
    d = D(c)*D(c)*D(c)
    # Set up a test array.
    test_array = [d, 0.1, 0.3, 0.5, 0.8, 0.999999, 1, 2, 7856432, 562345.983526627,
              987299826534888736635.981425, 99999999, '80996655413131.092562882',
              624262829789793071275019512095723496525655125151635435275759776234274,
              999999999999999999999999999999999999999999999999999999999999999999999,
              0.2345378907890453453456, 0.00000000000000000000000000001]
    # Set the methods to use.
    #method_array = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
    #                18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31]
    # List of well working algorithms.
    method_array = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                    13, 19, 20, 21, 23, 24, 25, 26, 31]
    # Perform some tests.
    for i in test_array:
        print("Test value:", i)
        for j in method_array:
            msgstr = "cubic_root_v" + str(j) + "():"
            print(msgstr)
            start = timer()
            prgstr = "cubic_root_v" + str(j) + "(D(" + str(i) + "))"
            cr = eval(prgstr)
            end = timer()
            print(cr)
            print("Elapsed time:", timedelta(seconds=end-start))
        print()
    # End of function. Return 1 for success.
    return 1

# Execute the script as module or as program.
if __name__ == '__main__':
    # Call the main script function.
    main()
