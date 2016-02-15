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
2. Define displacement steps sizes
3. Define a mean free function which gives l(r) 
4. Find the number of steps required to get out of shell given the current l(r).

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

speed_light = 3.0*10**10     #speed of light: unites: [cm / s] 

r_initial = 0
radius_fr = 0
scale_cuz_too_big = 10**0
R_final_cm   = scale_cuz_too_big * 7*10**10     #Solar radius in cm

R_final_fr = 1.0
print("The Initial Radius is 0.")

Mass_tot  = 2.0*10**33   #Mass of Sun Units: grams

trial_iterations = 1
num_steps = 5000         #Computes the number of steps to take in sun.  
dr_cm = R_final_cm/num_steps  #Defines the displacement scale.


time_list = []

r_0 = 0          #Sets initial radius to zero. 
radius_fr = 0 
x_pos = 0
y_pos = 0
z_pos = 0
distance_tot = 0
displacement_tot = 0
N_total = 0 

for i in range(0, num_steps):
	print("Start of loop. Iteration: %d" % i) 
	
	displacement_val = dr_cm
	MFP_val = get_MFP(radius_fr)
	N_val   = (displacement_val / MFP_val) ** 2.
	
	N_total += N_val	
	distance_tot += MFP_val * N_val      #Calculates the distance travelled in 
	

print("Total distance travelled: %.2e [cm]" % distance_tot)
time_val_sec = distance_tot / speed_light    #units: seconds
time_val_yrs = time_val_sec / 31622400.  #Convert seconds to years
	
print("It took this photon %.2e years to get out of the sun!" % time_val_yrs)
time_list.append(time_val_yrs)

		
time_total = 0
for k in range(0, len(time_list)):
	time_total += float(time_list[k])
		


time_average_yrs = time_total/len(time_list)

print("\n==================================================================== ")
print("*** The Time Average for a photon to exit the sun is %.2e years! *** " % time_average_yrs)


