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



print("Welcome this code will compute a quick and dirty random walk given a non-constant Mean Free Path.")

r_initial = 0
print("The Initial Radius is 0.")

Mass_tot  = 2.0*10**33   #Mass of Sun Units: grams
num_steps = 1000         #Computes the number of steps to take in sun.  

dr = Mass_tot/num_steps  #Defines the displacement scale.






