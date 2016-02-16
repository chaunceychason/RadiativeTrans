import sys
import os
import numpy as np
import math
import pandas as pd
from scipy import stats  


print("WELCOME TO A QUICK PYTHON PROGRAM.")
print("It will find roots given a single equation. yep. ")


"""
===============================
FUNCTION:
	
	(3-x)e^x  =  3  

 ===>  f(x) = (3-x)e^x - 3  =  0  

===============================
"""

def func( x ):
	#This is the function function its function is to be functionally functional
	return ((3. - x ) * math.exp(x) - 3.)


#Guess something
initial_guess = 1. 
delta = 0.0001
epsilon = 0.01

#Initialize relevant things
iterations = 0
max_iterations = 10000000


#Exclude the solution x=0 since it is trivial. 
exclude_value = 0

#Initialize Derivative
deriv = 0
initial_error = abs(func(initial_guess))

#First comput the guess and initial derivative.
if initial_error > epsilon: 
	print("First Guess return value:  %f" % initial_error)
	#Compute the derivative approximated using simple algorithm.
	deriv =  ((func(initial_guess+(delta/2.))) + func(initial_guess-(delta/2.))) / delta
	#Choose New Guess for x.
	change_value = 0.5 * initial_error / deriv 

#Compute f(x) after new guess.
new_guess = initial_guess -  ( change_value  )
while (abs(func(new_guess)) > epsilon) or (iterations>=max_iterations):
	#Ensures X doesnt get stuck in local min at zero. 
	if ((new_guess == 0) or (new_guess < -2.5) or (-0.0000001 < new_guess < 0.00000001)):
		if random.ranom() > 0.6:
			new_guess = -100. * random.random 
		else:
			new_guess = 100 * random.random

	error_value = func(new_guess)
	print("Guess: %.5f Error Value: %.4f" % (new_guess, error_value))
	deriv =  ((func(new_guess+(delta/2.))) + func(new_guess-(delta/2.))) / delta
	change_value = -error_value / deriv
	prev_guess = new_guess 
	new_guess += change_value
	iterations += 1


print("Exited the while loop. Either convergence, or max_iterations reached. ")
if max_iterations <= iterations:
	print("Bummer. Max Iterations Reached.")
else:
	print("YAY!. Printing values....\n\n")
	print("Roots for x:  {0.000 , %.3f } " % new_guess)


print("END. ")




#Doop. 