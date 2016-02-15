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
scale_cuz_too_big = 10**-8

R_final_cm   = scale_cuz_too_big * 7*10**10     #Solar radius in cm

R_final_fr = 1.0
print("The Initial Radius is 0.")

Mass_tot  = 2.0*10**33   #Mass of Sun Units: grams

trial_iterations = 1
num_steps = 2000         #Computes the number of steps to take in sun.  
dr_cm = 0.5*R_final_cm/num_steps  #Defines the displacement scale. units: [cm]


time_list = []

radius_fr = 0 
distance_tot = 0
displacement_tot = 0
N_total = 0 

#Compute the 1/2 radius to total radius:
displacement2 = 0.5*R_final_cm 
MFP2 = 1   #units: [cm]
N2 = (displacement2/MFP2)**2
distance_tot2 = N2 * MFP2 
time2_sec = distance_tot2/ speed_light
time2_years = time2_sec/31622400.


for i in range(0, num_steps):
	print("Start of loop. Iteration: %d" % i) 
	radius_fr += dr_cm/R_final_cm   #Flag: changed. 
	displacement_val = dr_cm
	MFP_val = get_MFP(radius_fr)
	N_val   = (displacement_val / MFP_val) ** 2.
	
	N_total += N_val	
	distance_tot += MFP_val * N_val      #Calculates the distance travelled in 
	

distance_total = distance_tot + distance_tot2
print("Total distance travelled: %.2e [cm]" % distance_total)
time_val_sec = distance_tot / speed_light    #units: seconds
time_val_yrs = time_val_sec / 31622400.  #Convert seconds to years

total_time = time_val_yrs + time2_years	
print("It took this photon %.2e years to get out of the sun!" % total_time)
time_list.append(total_time)


"""		
time_total = 0
for k in range(0, len(time_list)):
	time_total += float(time_list[k])
		
time_average_yrs = time_total/len(time_list)

print("\n==================================================================== ")
print("*** The Time Average for a photon to exit the sun is %.2e years! *** " % time_average_yrs)
"""

