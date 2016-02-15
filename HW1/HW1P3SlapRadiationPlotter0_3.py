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

1. Begin with j, and alpha, T, electron number density, wavelength. 
"""
	
def plot_Ivs(x_arr, y_arr):
	fignamestring = "IntensityvSlapJan25.png"

	fig = plt.figure(figsize=(6,6))
	#ax = fig.add_subplot(111)
	#ax.set_title("Intensity vs. Slab Distance" ,fontsize=14)
	#ax.set_xlabel("Slab Distance Traversed [cm]",fontsize=12)
	#ax.set_ylabel("Specific Intensity [cgs]",fontsize=12)
	#ax.grid(True,linestyle='-',color='0.75')

	#sc = plt.scatter(x_arr, y_arr, c=z_arr, s=sizearr, edgecolor='') 	
	for i in range(0,len(x_arr)):
		plt.scatter(float(x_arr[i]), float(y_arr[i])) 	
	
	#print x_arr
	plt.yscale('log')	
		

	plt.savefig(fignamestring)
	print("Saving Plot: %s " % fignamestring)
	plt.clf()
	return 0

def plot_xyz_displacement(displacement_arr, radius_arr):
	color_arr = np.linspace(0, max(displacement_arr), len(displacement_arr))
	plt.plot(color_arr, displacement_arr) 
	plt.title("Displacement Vs. Step")
	plt.xlabel("Step")
	plt.ylabel("Displacement [cm]")
	fignamestring = "DisplacementVstepAJan25.png"
	plt.savefig(fignamestring)
        print("Saving Plot: %s " % fignamestring)
	return 0

def conv_pc_to_cm(x_pc):
	x_cm = x_pc * 3.086*10**18.
	return x_cm

def calc_jff(n_e, T):
	jff = 1.4*10**-27. * (n_e*n_e) * (T**(-1./2.))
	jff = 0  
	print("WARNING: jff temporarily set to zero.")
	return jff

def calc_alphaff(n_e, T, freq):
	alphaff = 0.018 * (T**(-3./2)) * (n_e*n_e) * (freq ** -2.)
	return alphaff

print("Welcome this code will compute a quick and dirty random walk through a slab.")


speed_light = 3.0*10**10     #speed of light: unites: [cm / s] 
T_K = 10**4 
s = 0
s_final_pc = 10.
s_final_cm = conv_pc_to_cm(s_final_pc)

I_0_jy = 10  #Jy/steradian
I_0_cgs = I_0_jy * 10**-23. 

n_e = 1.0 #cm^-3

wavelength = 21.1  #cm
frequency = speed_light / wavelength 

num_steps = 200
ds = s_final_cm/num_steps 

s_array = []
I_array = []

while s < s_final_cm:
	if I_0_cgs < 0:
		print("YIKES! Intensity < 0. Error.")
		break

	a_ff = calc_alphaff(n_e, T_K, frequency)
	jff  = calc_jff(n_e, T_K) 
	dI_ds = -a_ff* I_0_cgs + jff 	
	I_1 = I_0_cgs + dI_ds * ds
	
	s_array.append(s)
	I_array.append(I_1)	

	I_0_cgs = I_1
	s += ds	

print("EXITED THE WHILE LOOP!")


print I_array
plot_Ivs(s_array, I_array)

#plot_xyz_displacement(displacement_arr, radius_arr)



