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


	
def plot_xyz_positions(x_arr, y_arr, z_arr):
	fignamestring = "X_YpositionsJan25.png"
	steparr = np.linspace(0, len(x_arr),len(x_arr))

	fig = plt.figure(figsize=(10,6))
	ax = fig.add_subplot(111)
	ax.set_title("X_Y Positions" ,fontsize=14)
	ax.set_xlabel("X",fontsize=12)
	ax.set_ylabel("Y",fontsize=12)
	ax.grid(True,linestyle='-',color='0.75')

	min_zarr = min(z_arr)
        max_zarr = max(z_arr)
	sizearr = []
	for i in range(0, len(z_arr)):
		sizeval = z_arr[i]+min_zarr
		sizeval = (sizeval / (max_zarr+min_zarr)) + 1.0
		sizearr.append(2.**sizeval)	
	
	
	#sc = plt.scatter(x_arr, y_arr, c=z_arr, s=sizearr, edgecolor='') 	
	sc = plt.scatter(x_arr, y_arr, c=steparr, s=sizearr, edgecolor='') 	
	#print x_arr

	plt.colorbar(sc)

	#ax.x_lim(min(x_arr),max(x_arr))
	#ax.y_lim(min(y_arr),max(y_arr))
	
	#plt.xlim(min(x_arr),max(x_arr))
	#plt.ylim(min(y_arr),max(y_arr))
	#plt.colorbar(z_arr)
	
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
	fignamestring = "DisplacementVstepJan25.png"
	plt.savefig(fignamestring)
        print("Saving Plot: %s " % fignamestring)
	return 0


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
#num_steps = 1000         #Computes the number of steps to take in sun.  
#dr = Mass_tot/num_steps  #Defines the displacement scale.

time_list = []
x_arr = []
y_arr = []
z_arr = []
displacement_arr = []
radius_arr = []

for i in range(0, trial_iterations):
	r_0 = 0          #Sets initial radius to zero. 
	radius_fr = 0 
	x_pos = 0
	y_pos = 0
	z_pos = 0
	distance_tot = 0
	iteration = 1
	print("Start of loop. Iteration: %d" % i) 
	p001radfr = True	
	p01radfr  = True
	p1radfr   = True
	p2radfr   = True
	p5radfr   = True
	p75radfr  = True
	p9radfr   = True

	while radius_fr < R_final_fr:
		random_xyz_array = np.random.rand(3,1)
		initial_random_xyz_array = random_xyz_array   #Values from 0 to 1. 

		random_xyz_array = random_xyz_array * 2.0 - 1.0    #operation to get values from -1 to 1. 
		#random_xyz_array * (3.**(1./2)) #Operation to get values from -1/3 to 1/3 so vector length is 1. 		

		x_val_init = random_xyz_array[0]    #values between -1 and 1.  Thus the max vector length is 3. 
		y_val_init = random_xyz_array[1]
		z_val_init = random_xyz_array[2]

		MFP_val = get_MFP(radius_fr)	#Get the Mean Free Path for this radius. 
		distance_init = x_val_init**2. + y_val_init**2. + z_val_init**2.  #Calculate the vector length for values.
		
		#Normalize the vector length to the MFP_length essentially randomizing direction. 
		alpha_norm = (MFP_val**2./3.0)**(1./2.)  #Normalizes the vector length
		alpha_norm * random_xyz_array            #Updates the xyz array with normalization.
		
		if iteration == 1:
			print initial_random_xyz_array
			print random_xyz_array
		
		x_val = float(random_xyz_array[0])    #values between -MFP and MFP.  Thus the max vector length is MFP
                y_val = float(random_xyz_array[1])
                z_val = float(random_xyz_array[2])

		x_pos += x_val
		y_pos += y_val
		z_pos += z_val
		
		x_arr.append(x_pos)
		y_arr.append(y_pos)
		z_arr.append(z_pos)
		
		distance = (x_val**2. + y_val**2. + z_val**2.)**(1./2.)	#Units: cm
		displacement = (x_pos**2 + y_pos**2 + z_pos**2)**(1./2.)  	#units: cm
		radius_fr = displacement/R_final_cm 		#units: solar radii 
		
		displacement_arr.append(displacement)
		radius_arr.append(radius_fr)	
		#print displacement
		
		if radius_fr > 0.001 and p001radfr == True:
			p001radfr = False
			print("Radius fraction is above 0.001.")
		if radius_fr > 0.01 and p01radfr == True:
			p01radfr = False
                        print("Radius fraction is above	0.01.")
		if radius_fr > 0.1 and p1radfr == True:
			p1radfr = False
                        print("Radius fraction is above 0.1.")
		if radius_fr > 0.2 and p2radfr == True:
			p2radfr = False
                        print("Radius fraction is above 0.2.")
		if radius_fr > 0.5 and p5radfr == True:
			p5radfr = False
                        print("Radius fraction is above 0.5.")
		if radius_fr > 0.75 and p75radfr == True:
			p75radfr = False
                        print("Radius fraction is above 0.75.")
		if radius_fr > 0.9 and p9radfr == True:
			p9radfr = False
                        print("Radius fraction is above 0.9.")

		
		
		distance_tot += distance    #add the travelled distance to the total distance. 							
		iteration += 1
	
	print("Finished the random walk out of the sun! Total Iterations: %d " % iteration)
	print("Total distance travelled: %.2f [cm]" % distance_tot)
	time_val_sec = distance_tot / speed_light    #units: seconds
	time_val_yrs = time_val_sec / 31622400.  #Convert seconds to years
	
	print("It took this photon %.2e years to get out of the sun!" % time_val_yrs)
	time_list.append(time_val_yrs)
	
	
plot_xyz_positions(x_arr, y_arr, z_arr)

plot_xyz_displacement(displacement_arr, radius_arr)

time_total = 0
for k in range(0, len(time_list)):
	time_total += float(time_list[k])
		

time_average_yrs = time_total/len(time_list)

print("\n==================================================================== ")
print("*** The Time Average for a photon to exit the sun is %.2e years! *** " % time_average_yrs)


