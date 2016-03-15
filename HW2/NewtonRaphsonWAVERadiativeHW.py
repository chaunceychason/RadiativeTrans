import sys
import os
import numpy as np
import math
import random
from scipy import stats  
import matplotlib
matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt
from datetime import datetime
startTime = datetime.now()

print("WELCOME TO A QUICK PYTHON PROGRAM.")
print("It will find roots given a single equation. yep. ")


"""
===============================
FUNCTION:
	
	(5-x)e^x  =  5  

 ===>  f(x) = (5-x)e^x - 5  =  0  

===============================
"""

def func( x ):
	#This is the function function its function is to be functionally functional
	return ((5. - x ) * math.exp(x) - 5.)

def calc_wavelength( x_val ):
	#Finds CGS wavelength given a root for x. 
	h = 6.62607 * (10 ** -27. )
	c = 3.0 * ( 10 ** 10.)
	k = 1.380648 * (10 ** -16.)
	T = 533.15

	wave = h * c / ( k * T * x_val )
	return wave


def plot_basic( xlist, ylist, title, xlab, ylab, legend_val, psize, xlog, ylog, yflip , pcounter):
	print("Entered Basic Plot Function")
	
	if len(xlist) != len(ylist):
		print("ERROR! X and Y DATA LENGTHS ARE DIFFERENT!")
		print("Length: x_data: %g" % len(xlist))
		print("Length: y_data: %g" % len(ylist))
	if len(xlist)==0 or len(ylist)==0:
		print("ERROR: list length is ZERO!")
		print("Length: x_data: %g" % len(xlist))
		print("Length: y_data: %g" % len(ylist))
	else:
		print("Length: x_data: %g" % len(xlist))
		print("Length: y_data: %g" % len(ylist))

	if legend_val != 0:
		pass

	plot_title=" Blank Title "   
	x_axis="Blank X"
	y_axis="Blank Y"
	pointsize = 5
	#figure_name=os.path.expanduser('~/Feb9LSSHW2_plot1_A' +'.png')
	#figure_name=os.path.expanduser('~/Feb9LSSHW2_plot1_A' +'.png')
	#Choose which type of plot you would like: Commented out.

	#sets new plot features from call. 	
	if True:
		plot_title = title
		x_axis = xlab
		y_axis = ylab
 		pointsize = psize

	
	plt.scatter(xlist, ylist, s=pointsize, lw=0)
	plt.title(plot_title)
	plt.xlabel(x_axis)
	plt.ylabel(y_axis)
	#plt.yscale("log")
	
	if yflip == True:
		try:
			plt.ylim(max(ylist), min(ylist))
		except:
			print("uh.oh.... try except statement. check ylim.")

	if ylog != 0:
		plt.yscale("log", nonposy='clip')
	if xlog != 0:
		plt.xscale("log", nonposy='clip')

	#plt.xlim(min(xlist), max(xlist) )
	
	print	
	figure_name=os.path.expanduser('~/Feb16RadTransHW2NR_plot%s.png' % pcount)

	plt.savefig(figure_name)
	print("Saving plot: %s" % figure_name)
	plt.clf()

	dummy = pcounter + 1
	return dummy



#Guess something
#initial_guess = 1
initial_guess = -1. 
delta = 0.001
epsilon = 0.001

#Initialize relevant things
iterations = 0
max_iterations = 10000000
guessesLIST = []
errorsLIST  = []

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
	change_value = 0.75 * initial_error / deriv 

guessesLIST.append(initial_guess)
errorsLIST.append(initial_error)

#Compute f(x) after new guess.
new_guess = initial_guess -  ( change_value  )
while (abs(func(new_guess)) > epsilon) or (iterations>=max_iterations):
	#Ensures X doesnt get stuck in local min at zero. 
	if ((new_guess == 0) or (new_guess < -2.5) or (-0.01 < new_guess < 0.01)):
		if random.random() > 0.6:
			new_guess = -20. * random.random() 
		else:
			new_guess = 20 * random.random()

	error_value = func(new_guess)
	print("Guess: %.5f Error Value: %.4f" % (new_guess, error_value))
	deriv =  ((func(new_guess+(delta/2.))) + func(new_guess-(delta/2.))) / delta
	change_value = -error_value / deriv
	prev_guess = new_guess 
	
	if (abs(error_value) < (10 * epsilon)):
		new_guess += (0.25 * change_value)
	else:
		new_guess += change_value

	guessesLIST.append(prev_guess)
	errorsLIST.append(error_value)
	
	iterations += 1


print("Exited the while loop. Either convergence, or max_iterations reached. ")
if max_iterations <= iterations:
	print("Bummer. Max Iterations Reached.")
else:
	print("YAY!. Printing values....\n\n")
	print("Roots for x:  {0.0000 , %.4f } " % new_guess)

print("Calculating Wavelength Max...")

wavelength_value = calc_wavelength( new_guess )
print("Wavelength [cm]: %.5e " % wavelength_value)
print("Wavelength [A] : %.5e " % (wavelength_value*(10**8.)))


print("Now Plotting. ...")

"""
===================
PLOTTING GUESSES 
===================
"""
pcount = 0
title_label = "Function Value F(x) vs. Guess" 
x_label = "Iteration"
y_label = "F(x)"
x_data  = np.linspace(0, len(errorsLIST), num=len(errorsLIST))
y_data  = errorsLIST
legend_val = 0
pointsize = 3
yflip = False
ylog = 0
xlog = 0

pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, \
         legend_val, pointsize, xlog, ylog,  yflip, pcount)


title_label = "Guess vs. Iteration" 
x_label = "Iteration"
y_label = "Guess"
x_data  = np.linspace(0, len(errorsLIST), num=len(errorsLIST))
y_data  = guessesLIST
legend_val = 0

#COMPUTE ERROR_SCALELIST
error_scaleLIST = []
max_error = abs(max(errorsLIST))
davalue_int = 1.
for idx, val in enumerate(errorsLIST):
	davalue = (val / max_error) * 25.
	davlue_int = int(davalue)
	if davalue < 1.:
		davalue_int == 1.
	error_scaleLIST.append(davalue_int)

pointsize = error_scaleLIST
yflip = False
ylog = 0
xlog = 0

pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, \
         legend_val, pointsize, xlog, ylog,  yflip, pcount)



print("Time: ")
print (datetime.now() - startTime)

print("END. ")




#Doop. 
