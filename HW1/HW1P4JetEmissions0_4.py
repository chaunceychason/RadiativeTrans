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
	
def plot_jvtheta(x_arr, y_arr, beta):
	fignamestring = "JetEmissionsBeta_%.1eJan26.png" % beta 

	fig = plt.figure(figsize=(6,6))
	ax = fig.add_subplot(111)
	title_string = "Volume Emissivity v. Theta Beta: %.1e" % beta
	ax.set_title(title_string ,fontsize=14)
	ax.set_xlabel("Theta",fontsize=12)
	ax.set_ylabel("j [units of j0]",fontsize=12)
	ax.grid(True,linestyle='-',color='0.75')
	"""
	#sc = plt.scatter(x_arr, y_arr, c=z_arr, s=sizearr, edgecolor='') 	
	norm_x = []
	norm_y = []
	max_x = max(x_arr)
	min_x = min(x_arr)
	max_y = max(y_arr)
	min_y = min(y_arr)

	print max_x
	print min_x
	print max_y
	print min_y
	
	for i in range(0,len(x_arr)):
		norm_x.append((x_arr[i])/max_x)
		norm_y.append((y_arr[i])/max_y)
	
	"""	
	#plt.plot(norm_x, norm_y) 	
	plt.plot(x_arr, y_arr) 	
	
	#print x_arr
	plt.yscale('log')	
		

	plt.savefig(fignamestring)
	print("Saving Plot: %s " % fignamestring)
	plt.clf()
	return 0


print("Welcome this code will compute a quick and dirty random walk through a slab.")
speed_light = 3.0*10**10     #speed of light: unites: [cm / s] 


beta1 = 1*10**-4
beta2 = 1*10**-2
beta3 = 0.99 
beta4 = 1.0

beta = beta4

j_array = []
x_array = []

x_nparr = np.linspace(0,2*np.pi, 200)
	
for element in x_nparr:
	j_array.append(1./(1. - beta*math.cos(element))**4.)
	

plot_jvtheta(x_nparr, j_array, beta)



