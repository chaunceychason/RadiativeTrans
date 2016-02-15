import numpy as np
import matplotlib
matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt



x_array = [10**-8, 10**-7, 10**-6, 10**-5, 10**-4, 10**-3.5, 10**-3]
y_array = [9*10**-16 , 9*10**-15, 8*10**-14, 8*10**-13, 1*10**-11, 2.1*10**-11, 1.13*10**-10]
x2_array = [10**-8,        10**-7,       10**-5,     10**-3,      10**-2,     10**-1,       10**0]
y2_array = [1.3*10**-13 , 1.3*10**-11, 1.3*10**-7, 1.3*10**-3, 1.3*10**-1, 1.29*10**1, 1.29*10**3]


plt.plot(x_array, y_array, "bo", markersize = 10, linestyle='-', rasterized=True)
plt.plot(x2_array, y2_array, "go", markersize = 10, linestyle='-', rasterized=True)
#plt.ylim((0,130))
#plt.xlim((0,130))	

tmp_filename = "gloriousplotJan24.png"
x_axis = "Total Radius Fraction" 
y_axis = "Estimated Random Walk Time [years]"
plot_title = "Fractional Radius vs. Random Walk Time."

#Sets plot parameters
plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)
plt.yscale('log')	
plt.xscale('log')

#Saves Plot
plt.savefig(tmp_filename)

