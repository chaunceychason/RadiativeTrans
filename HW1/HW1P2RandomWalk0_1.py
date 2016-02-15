import sys
import os
import math
import numpy as np
#from scipy import interpolate
import matplotlib
matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt
from random import randint


"""
Code Design:

1. Begin with inital radius and Mass of star. 
2. Define displacement steps sizes. 
3. Define a mean free function which gives l(r) 
4. Choose a random direction to step in, and update current position.
	using: randint(0,360)

5. Record the total distance travelled and iterate to next step.
6. Stop once youved reached R total displacement. 

"""
def get_MFP(radius):
	#Using Linear Function
	if radius < 0.5:
		MFP_val = 2.0*radius+0.01	#units: cm
	else: 
		MFP_val = 1.0			#units: cm
	return MFP_val


print("Welcome this code will compute a quick and dirty random walk given a non-constant Mean Free Path.")

r_initial = 0
R_final_cm   = 7*10**10     #Solar radius in cm
R_final_fr = 1.0
print("The Initial Radius is 0.")

Mass_tot  = 2.0*10**33   #Mass of Sun Units: grams
trial_iterations = 1
num_steps = 1000         #Computes the number of steps to take in sun.  


dr = Mass_tot/num_steps  #Defines the displacement scale.

for i in range(0, trial_iterations):
	r_0 = 0          #Sets initial radius to zero. 
	x_pos = 0
	y_pos = 0
	z_pos = 0
	distance_tot = 0
	iteration = 1
	print("Start of loop. Iteration: %d" % i) 
	
	while radius < R_final_fr:
		random_xyz_array = np.random.rand(3,1)
		initial_random_xyz_array = random_xyz_array
		random_xyz_array * 2.0 - 1.0  #operation to get values from -1 to 1. 
		x_val_init = random_xyz_array[0]
		y_val_init = random_xyz_array[1]
		z_val_init = random_xyz_array[2]

		
		if iteration = 1:
			print initial_random_xyz_array
			print random_xyz_array
		
		x_pos += x_val
		y_pos += y_val
		z_pos += z_val

		distance = x_val**2. + y_val**2. + z_val**2.
				


		iteration += 1






