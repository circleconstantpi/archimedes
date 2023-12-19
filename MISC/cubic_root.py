#!/usr/bin/python3
# -*- coding: utf-8 -*-
'''Collection of algorithm for the calculation of the cubic root.

Following famous algorithms are implemented:

Regula Falsi Method, Method of False Position, False Position Method
Bisection Method
Binary Search Method
Newton's Method
Newton–Raphson Method
Halley's Method
Householder's Method
Steffensen Method**

** works only for small values as expected

Implementation of various procedures for study purposes.

USAGE: (from shell prompt)
    python3 ./cubic_root.py

See also:
www.math-cs.gordon.edu/courses/mat342/python/findroot.py
cocalc.com/share/public_paths/embed/bf8c2613e8e9f04f4b052a83fb30a14e1ccea697
'''
# pylint: disable=invalid-name
# pylint: disable=redefined-outer-name
# pylint: disable=no-else-return
# pylint: disable=chained-comparison
# pylint: disable=invalid-unary-operand-type
# pylint: disable=line-too-long
# pylint: disable=eval-used

__author__ = "Dr. Peter Netz"
__copyright__ = "Copyright (C) 2023, Dr. Peter Netz"
__license__ = "MIT"
__version__ = "0.4"

# Import the standard Python modules.
from timeit import default_timer as timer
from datetime import timedelta

# Import the standard Python modules.
from decimal import Decimal as D
from decimal import getcontext, localcontext, InvalidOperation, ROUND_HALF_EVEN, ROUND_HALF_DOWN

# Initialise the constant.
PRECISION = 128

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
    xn = (2*x0 + D(a)/D(x0*x0)) / 3
    # Change the local context.
    with localcontext() as ctx:
        # Change the local context behaviour.
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
        ctx.rounding = ROUND_HALF_DOWN
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
# Function cubic_root_v9()
# ----------------------------------------------------------------------
def cubic_root_v9(num):
    '''Applying Regula falsi to get the cubic root.'''
    # Definition of governing equation.
    def f(x, a):
        return x*x*x - a
    # Calculate the number of leading digits.
    cln = len(str(num).split(".")[0])
    # Get the used decimal precision.
    c = getcontext()
    prec = c.prec-cln
    # Set the convergence criterion.
    # eps: half of upper bound for relative error.
    eps = D(10)**(-prec)
    # Initialise the variable side.
    side = 0
    # a,b are endpoints of an interval where we are searching.
    a = -2 * num
    b = 2 * num
    # Starting values at the endpoints of the interval.
    fa = f(a, num)
    fb = f(b, num)
    # Perform a calculation in a new block using a local context.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding = ROUND_HALF_DOWN
        # Run an infite loop.
        while True:
            c = (fa * b - fb * a) / (fa - fb)
            if abs(b - a) < eps * abs(b + a):
                break
            fc = f(c, num)
            if fc * fb > 0:
                # fc and fb have the same sign, copy c to b.
                b = c
                fb = fc
                if side == -1:
                    fa /= 2
                side = -1
            elif fa * fc > 0:
                # fc and fa have the same sign, copy c to a.
                a = c
                fa = fc
                if side == +1:
                    fb /= 2
                side = +1
            else:
                # fc * fa, fb is very small (or looks like zero).
                break
    # Restore precision.
    c = +c
    # Return the cubic root.
    return c

# ----------------------------------------------------------------------
# Function cubic_root_v10()
# ----------------------------------------------------------------------
def cubic_root_v10(a):
    '''Steffensen method.'''
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
    p0 = D(1)/D(2)
    # Run calculation in a block.
    with localcontext() as ctx:
        # Change the local context.
        ctx.prec += 2
        ctx.rounding=ROUND_HALF_DOWN
        # Run an infinite loop.
        while True:
            # Try to calculate the cubic root.
            try:
                p1 = p0 + f(p0, a)
                p2 = p1 + f(p1, a)
                p = p2 - (p2 - p1)*(p2 - p1) / (p2 - 2*p1 + p0)
                # Convergence criteria.
                if abs(p-p0) < eps:
                    break
                # Store value for next loop.
                p0 = p
            except InvalidOperation:
                break
    p = +p
    return p

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

#method_array = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
method_array = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

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
