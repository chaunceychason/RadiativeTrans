import sys
import math
import numpy as np
import matplotlib
matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
def calc_gamma(v):
	c = 3. * 10**8.
	gamma = (1./(1. - (v**2/c**2))**0.5)
	return gamma

def calc_a(q, E, v, m):
	gamma = calc_gamma(v)
	a = ((q*E)-((q*v*v*E)/(c*c)))/(gamma*m) 
	return a


#THIS PROGRAM PLOTS DATA IN FORM 'x, y'
print("WELCOME. ")

plot_title="r(t) vs t "
x_axis="Time [sec] "
y_axis=" r ( t )  [meters]"

plt.title(plot_title)
plt.xlabel(x_axis)
plt.ylabel(y_axis)

end_seconds = 0.004
num_points  = 1000

time_array = [end_seconds * float(i)/float(num_points)  for i in range(0, num_points)]
delta_t    = time_array[1] - time_array[0]
value_array = [0 for i in range(0, num_points)]
beta_array  = [0 for i in range(0, num_points)]

c = 3.0e8   # m/s
q = 1.6e-19 # charge 
m = 9.1e-31 # kg
E_o = 1.0  
const1 = ( q * E_o / m )
E = E_o
beta = 0
beta_bool = 0

x_o = 0
v_o = 0
a_o = 0
x   = x_o
v   = v_o
a = calc_a(q, E, v, m)	
print a
print delta_t

for idx, time in enumerate( time_array ):
	value_array[idx] = x 
	v = v_o + a * delta_t
	x = (x_o + v*delta_t + 1./2.* a * delta_t*delta_t)
	a = calc_a(q, E, v, m)
	x_o = x 
	v_o = v 
	#value = (c * (  (c*c + const1*const1 * time * time )**(0.5) ) - c ) / ( const1 ) 
	#value_array[idx] = value
	beta = v / c 
	if beta > 0.5 and beta_bool == 0:
		beta_bool += 1
		beta_idx = idx
	print("Beta = %.4f" % beta)

c
plt.plot(time_array, value_array ,color='g' , label='Electron in Eo' )
plt.plot(time_array[beta_idx], value_array[beta_idx], color='r', marker='D', label='Beta = 0.5')


#---------------	
#plt.ylim(0, 1400)	
plt.legend(loc='upper left')

#Saves Plot
tmp_filename = "plottrajectoryMar15.png"
plt.savefig(tmp_filename, rasterized=True)
plt.clf()






print "The program has finished running. All files closed. \nThe results should be in your directory"



#End