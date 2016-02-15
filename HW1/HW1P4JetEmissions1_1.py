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
	fignamestring = "JetEmissionsBetasJan26.png" 

	fig = plt.figure(figsize=(6,6))
	ax = fig.add_subplot(111)
	#title_string = "Volume Emissivity v. Theta by Beta"
	#ax.set_title(title_string ,fontsize=14)
	ax.set_xlabel("Theta",fontsize=12)
	ax.set_ylabel("j [units of j0]",fontsize=12)
	ax.grid(True,linestyle='-',color='0.75')

	plt.plot(x_arr, y_arr, label=str(beta)) 	
	plt.yscale('log')	
		
	#plt.savefig(fignamestring)
	#print("Saving Plot: %s " % fignamestring)
	#plt.clf()
	return 0

print("Welcome this code will compute a quick and dirty random walk through a slab.")
speed_light = 3.0*10**10     #speed of light: unites: [cm / s] 

beta1 = 1*10**-4
beta2 = 1*10**-2
beta3 = 1*10**-1
beta4 = 0.5
beta5 = 0.75
beta6 = 0.9
beta7 = 1.0 

beta_array = [beta1, beta2, beta3, beta4, beta5, beta6, beta7]
#beta = beta4

x_nparr = np.linspace(0,2*np.pi, 200)


for i in range(0, len(beta_array)):

	beta = float(beta_array[i])
	j_array = []
	x_array = []

	i = 0	
	for element in x_nparr:	
		if (1.-beta*math.cos(element))!=0:
			j_array.append(1./(1. - beta*math.cos(element))**4.)
			x_array.append(x_nparr[i])
			i += 1
		else:
			i += 1
			continue

	plot_jvtheta(x_array, j_array, beta)

	del j_array[:]
	del x_array[:]


plt.legend('best')

figurename = "MultiJetEmissionsJan26.png"
plt.savefig(figurename)
plt.clf()
