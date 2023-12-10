#!/usr/bin/python3

# Import the standard Python module math.
import math

# Define the startvalues.
a0 = 2*math.sqrt(3)
b0 = 3

# Set the number of iterations.
iteration = 4

# Loop an iteration from 0 to 5 to get 5 values of Pi.
for i in range(0, iteration+1):
    # In the first loop the startvalues gives Pi.
    if i == 0:
        newpi = (a0 + b0) /2
        print(newpi)
    else:
        # Calculate iterative Pi.
        a1 = (2*a0*b0) / (a0 + b0)
        b1 = math.sqrt(b0*a1)
        a0 = a1
        b0 = b1
        b3 = (3*a1*b1) / (2*a1+b1)
        a3 = (a1*b1**2)**(1/3)
        newpi = (a3 + b3) / 2
        print(newpi)
