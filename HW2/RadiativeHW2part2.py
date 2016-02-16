import sys
import os
import numpy as np
import math
import matplotlib
matplotlib.use( 'Agg' )
import matplotlib.pyplot as plt
from datetime import datetime
startTime = datetime.now()

def plot_basic( xlist, ylist, title, xlab, ylab, legend_val, psize, xlog, ylog, yflip , pcounter):
	print("Entered Basic Plot Function")
	
	if len(xlist) != len(ylist):
		print("ERROR! X and Y DATA LENGTHS ARE DIFFERENT!")
		print("Length: x_data: %g" % len(xlist))
		print("Length: y_data: %g" % len(ylist))
	if len(xlist)==0 or len(ylist)==0:
		print("ERROR: list length is ZERO!")
		print("Length: x_data: %g" % len(xlist))
		print("Length: y_data: %g" % len(ylist))
	else:
		print("Length: x_data: %g" % len(xlist))
		print("Length: y_data: %g" % len(ylist))

	if legend_val != 0:
		pass

	plot_title=" Blank Title "   
	x_axis="Blank X"
	y_axis="Blank Y"
	pointsize = 5
	#figure_name=os.path.expanduser('~/Feb9LSSHW2_plot1_A' +'.png')
	#figure_name=os.path.expanduser('~/Feb9LSSHW2_plot1_A' +'.png')
	#Choose which type of plot you would like: Commented out.

	#sets new plot features from call. 	
	if True:
		plot_title = title
		x_axis = xlab
		y_axis = ylab
 		pointsize = psize

	#plt.plot(xlist, ylist)
	#plt.plot(xlist, ylist, linestyle="", marker="o", markeredgecolor=None, markeredgewidth=0.0)
	#plt.plot(scale1, Vmaxarray, linestyle="-", marker="o")
	plt.scatter(xlist, ylist, s=pointsize, lw=0)
	plt.title(plot_title)
	plt.xlabel(x_axis)
	plt.ylabel(y_axis)
	#plt.yscale("log")
	
	if yflip == True:
		try:
			plt.ylim(max(ylist), min(ylist))
		except:
			print("uh.oh.... try except statement. check ylim.")

	if ylog != 0:
		plt.yscale("log", nonposy='clip')
	if xlog != 0:
		plt.xscale("log", nonposy='clip')

	plt.xlim(min(xlist), max(xlist) )
	#plt.ylim(min(ylist), max(ylist) )

	#plt.savefig(figure_name)
	#print("Saving plot: %s" % figure_name)
	print	
	figure_name=os.path.expanduser('~/Feb14RadTransHW2_plot%s.png' % pcount)

	plt.savefig(figure_name)
	print("Saving plot: %s" % figure_name)
	plt.clf()

	dummy = pcounter + 1
	return dummy


def plot_hist( xlist, num_of_bins, hist_weights, title, xlab, ylab, log_y, pcounter):
	print("Entered HISTOGRAM Plot Function")
	
	plot_title=" Blank Title "   
	x_axis="Blank X"
	y_axis="Blank Y"
	

	#print xlist
	print("==========================")
	print("type " + str(type(xlist)))
	print("len " + str(len(xlist)))

	if True:
		plot_title = title
		x_axis = xlab
		y_axis = ylab
 		

	plt.hist(xlist, bins=num_of_bins, weights=hist_weights)
		
	plt.title(plot_title)
	plt.xlabel(x_axis)
	plt.ylabel(y_axis)
	#plt.yscale("log")
	if log_y != 0:
		plt.yscale('log', nonposy='clip')

	print
	
	figure_name=os.path.expanduser('~/Feb14RadTransHW2_plot%s.png' % pcount)
	plt.savefig(figure_name)
	print("Saving plot: %s" % figure_name)
	plt.clf()
	#Tracks plot number to assign unique name. 
	dummy = pcounter + 1
	return dummy

def integrate_basic( BLIST, boundsLIST, low_bound, high_bound ):
	#Integrats the function numerically. 
	integral_value = 0
	"""
	#TOTAL BOUND INTEGRATION
	for i, boundval in enumerate(boundsLIST[:-1]):
		binwidth = abs(boundsLIST[i+1]-boundsLIST[i])
		current_value = BLIST[i] * boundval * (binwidth)
		integral_value += current_value
	return integral_value
	"""
	#PARTIAL BOUND INTEGRATION
	for i, boundval in enumerate(boundsLIST[:-1]):
		if (boundval >= (10.**low_bound)) and (boundval <= (10.**high_bound)):
			binwidth = abs(boundsLIST[i+1]-boundsLIST[i])
			current_value = BLIST[i] * boundval * (binwidth)
			integral_value += current_value
	return integral_value

def find_z_max_given_M( Mag_list, z_list,  M_r ):
	z_Max = 0 
	#Enumerate the redshifts. 
	for i, x in enumerate(z_LIST):
		#Find the Max redshift for a given M_r
		if Mag_list[i] == M_r and z_Max < x:
			z_Max = x  
	return z_Max

def get_volume(max_z):
	#This volume uses max_redshift (max_z) to compute volume of SDSS sample: 
	#	Note: Uses max_z, instead of median since it is a volume limited, instead of flux-lim. 
	"""
	#Note: 
		Could compute the volume accurately using the following resources:
  				http://home.fnal.gov/~gnedin/cc/    
		This site gives the distance between two redshifts. 

	#Alternative: Use the approximation that:  D ~= 3000 * z
	 This expression is actually more accurately D ~= 2950 * z. 	
	"""
	radius = 2950.*max_z 
	DR7_SAMPLE_COVERAGE_STERADIANS = 2.295   #Should be set to 2.295 sterad. (7675.2 deg^2)
	fract_vol = DR7_SAMPLE_COVERAGE_STERADIANS / (4.*np.pi)    #Fractional sky coverage DR7 to 4pi. 
	vol = fract_vol*(4./3.) * np.pi * (radius**3.) 
	return vol  #Units: [h-1 Mpc^3]

def get_surfacearea_mks( radius  ):
	#Calculates the volume given a radius. 
	return (4. * np.pi * (radius ** 2.))

def get_surfacearea_cgs( radius  ):
	#Calculates the volume given a radius. 
	return (4. * np.pi * (radius ** 2.))

def calc_B_freq( temp, freq):
	#Calculates the Specific Intensity given freq, Units: W sr-1 m-2 Hz-1
	"""
	h_cgs = 6.6261 * (10**-27.)
	c_cgs = 3.0 * (10 ** 10. )
	k_bolz = 1.3806 * ( 10 ** -16. )
	"""
	h_mks = 6.626 * (10**-34.)
	c_mks = 3. * (10**8.)
	k_bolz = 1.38 * (10**-23.)

	B_num   = (2. * h_mks * (freq ** 3.))
	try:
		B_denom = (c_mks ** 2.) * ( math.exp(h_mks*freq/(k_bolz*temp)) - 1.0)
		B_freq =  B_num / B_denom
	except:
		B_freq = 0.
	#if B_freq < (10**-18):
	#	B_freq = 0.
	return B_freq

def calc_B_wavelength( temp, wave):
	#Calculates the Specific Intensity given wavelength, Units: W sr-1 m-3
	#h_ergs = 1.05457 * (10**-34.)  
	#h_cgs = h_ergs * (10**7.)
	"""
	h_cgs = 6.6261 * (10**-27.)
	c_cgs = 3.0 * (10 ** 10. )
	k_bolz = 1.3806 * ( 10 ** -16. )
	"""
	h_mks = 6.626 * (10**-34.)
	c_mks = 3. * (10**8.)
	k_bolz = 1.38 * (10**-23.)

	B_num   = (2. * h_mks * (c_mks ** 2. ))
	try:
		B_denom = ( wave ** 5.) * ( math.exp(h_mks*c_mks/(wave*k_bolz*temp)) - 1.0)
		B_wave =  B_num/ B_denom
	except:
		B_wave = 0.
	#if B_wave < (10**-4):
	#	B_wave = 0.
	return B_wave 

def calc_B_freq_cgs( temp, freq):
	#Calculates the Specific Intensity given freq, Units: W sr-1 m-2 Hz-1
	
	h_cgs = 6.6261 * (10**-27.)
	c_cgs = 3.0 * (10 ** 10. )
	k_bolz = 1.3806 * ( 10 ** -16. )
	"""
	h_mks = 6.626 * (10**-34.)
	c_mks = 3. * (10**8.)
	k_bolz = 1.38 * (10**-23.)
	"""
	B_num   = (2. * h_cgs * (freq ** 3.))
	try:
		B_denom = (c_cgs ** 2.) * ( math.exp(h_cgs*freq/(k_bolz*temp)) - 1.0)
		B_freq =  B_num / B_denom
	except:
		B_freq = 0.
	#if B_freq < (10**-18):
	#	B_freq = 0.
	return B_freq

def calc_B_wavelength_cgs( temp, wave):
	#Calculates the Specific Intensity given wavelength, Units: W sr-1 m-3
	#h_ergs = 1.05457 * (10**-34.)  
	#h_cgs = h_ergs * (10**7.)
	
	h_cgs = 6.6261 * (10**-27.)
	c_cgs = 3.0 * (10 ** 10. )
	#k_bolz = 1.3806 * ( 10 ** -16. )
	#ergs cm-2 K-4 s-1
	k_bolz = 1.3806 * ( 10 ** -16. ) 
	"""
	h_mks = 6.626 * (10**-34.)
	c_mks = 3. * (10**8.)
	k_bolz = 1.38 * (10**-23.)
	"""
	B_num   = (2. * h_cgs * (c_cgs ** 2. ))
	try:
		B_denom = ( wave ** 5.) * ( math.exp(h_cgs*c_cgs/(wave*k_bolz*temp)) - 1.0)
		B_wave =  B_num/ B_denom
	except:
		B_wave = 0.
	#if B_wave < (10**-12):
	#	B_wave = 0.
	return B_wave 
def calc_weins_wave( T_k ):
	#Returns back in meters given T in kelvin; b [units; hertz/ kelvin]
	b_constwave = 2.89777 * (10**-3) 
	wave_max_weins_meters = b_constwave / (T_k)
	wave_max_weins_cm = wave_max_weins_meters * 100.
	return wave_max_weins_cm

def calc_weins_freq( T_k ):
	#Returns back in meters given T in kelvin; b [units; hertz/kelvin]
	b_constfreq = 5.8789 * (10.**10.)
	freq_max_weins = b_constfreq * (T_k)
	return freq_max_weins


def	calc_wavelength_given_energy(energy_eV ):
	#Returns Value in cm. 
	h_const = 6.626 * (10**-34.)
	energy_joule = energy_eV * (1.6022*(10**-19.))
	c_mks = 3.0 * (10**8.)
	wavelength_value_mks = h_const * c_mks / (energy_joule)
	wavelength_value_cgs = wavelength_value_mks * 100.
	return wavelength_value_cgs


def edd_luminosity( mass_sol ):
	#Returns Units in Watts. may be called otherwise. 
	edd_luminosity.solarlumin = 3.20 * (10.**4. ) * mass_sol

	edd_luminosity.watts      = 1.26 * (10.**31.) * mass_sol 
	return edd_luminosity.watts

def conv_eV_cm( eV_value ):
	return 

def conv_F_to_K( t_F ):
	return float((t_F + 459.67) * (5./9.))

def conv_watt_to_ergPERs( watt_value):
	return watt_value * 10000000

def checkcondition(condition1, condition2):
	if (condition1 or condition2) == False:
		print("Conditions evaluated to FALSE! This should break. ")
		return False
	else:
		#The conditions are true!
		return True

"""
================================================================
================================================================
	 

   
HOMEWORK #2.  Written by:    Nicholas Chason 
----------------------------------------------------------------
INSTRUCTIONS: Integrates and plots the Blackbody spectrum. 
A.  - i.  
	- ii. 

================================================================
================================================================
"""
#CONSTANTS:
h_cgs = 6.6261 * (10**-27.)
c_cgs = 3.0 * (10 ** 10. )


print("Number of arguments: %s" % len(sys.argv))
print("Argument List:", str(sys.argv))
#boundLIST  = []
B_wave_mksLIST = []
B_freq_mksLIST = []
B_wave_cgsLIST = []
B_freq_cgsLIST = []
B_waveNORMLIST = []
B_freqNORMLIST = []
#Kelvin
Temperture_F = 100. 
#Bins
num_bins = 1000
"""
=====================
ENERGY SETTINGS
=====================
"""

#Energy in eV 
specific_energy = 2500.  
#Corresponds to 5000 angstroms wavelength , energy 
#specific_energy = 2.479  

#In CM
specific_wave = calc_wavelength_given_energy(specific_energy)

if len(sys.argv) == 1:
	print("setting default values.")
if len(sys.argv) == 2:
	#Kelvin
	Temperture_F = float(sys.argv[1])
	print("Temperture: %.2f F" % (Temperture_F))
"""
if len(argv) == 3:
	Temperture = sys.argv[1] 
	print("Temperture: %.2f " % Temperture)
"""
#Records key values: Peaks are in linear space. 
Temp_k = conv_F_to_K( Temperture_F )
wave_peak = calc_weins_wave( Temp_k )
freq_peak = calc_weins_freq( Temp_k )

#Set Bounds in cm; in log space 
#lower_bound_wave = -9.
#upper_bound_wave = 1.
wave_peak_log = np.log10(wave_peak)
lower_bound_wave =wave_peak_log -2.
upper_bound_wave = wave_peak_log+2.


#Set Bounds in cm 
#lower_bound_freq = 9.
#upper_bound_freq = 16.
freq_peak_log =  np.log10(freq_peak)
lower_bound_freq = freq_peak_log-2.
upper_bound_freq = freq_peak_log+2.

print("\nPEAKS:\n Wave: %.2f \n Freq: %.2f \n" % (wave_peak_log, freq_peak_log))
print("lower wave: %.4f" % lower_bound_wave)
print("upper wave: %.4f" % upper_bound_wave)
print("lower freq: %.4f" % lower_bound_freq)
print("upper freq: %.4f \n" % upper_bound_freq)


print("SETTING EYE VALUES!")
specific_energy_eye = wave_peak
specific_wave_eye = calc_wavelength_given_energy(specific_energy_eye)


boundwaveLIST = np.logspace(lower_bound_wave, upper_bound_wave, num_bins)
boundfreqLIST = np.logspace(lower_bound_freq, upper_bound_freq, num_bins)



"""
#-------------------------------------------------
# OPENING FILE: Sets the datafile name for opening
#-------------------------------------------------
datafilename = "./SDSS_DR7ordered.dat"     #An ordered version to save computation. 
ordered_file = True


if os.path.exists(datafilename):
	print("\n\t\tCONGRATS!\nSuccessfully opened %s ..." % datafilename)
	pass
else:
	print("os.path.exists says: FILE DOES NOT EXIST....")
	print("QUITTING PROGRAM.")
	sys.exit("Error Encountered with File I/O")    
"""
#plot number counter
pcount = 1

#Variable that stores the index of the subsample cutoff point
current_index = 0

"""
=====================================================
Calculate Intensity Lists for given Bounds
=====================================================
"""
#In MKS Units, then CGS Units
for idx, wave_value in enumerate(boundwaveLIST):
	B_wave = calc_B_wavelength( Temp_k, wave_value )
	B_wave_mksLIST.append(B_wave)
	B_wave_cgs = calc_B_wavelength_cgs( Temp_k, wave_value )
	B_wave_cgsLIST.append(B_wave_cgs)

for idx2, freq_value in enumerate(boundfreqLIST):
	B_freq = calc_B_freq( Temp_k, freq_value )
	B_freq_mksLIST.append(B_freq)
	B_freq_cgs = calc_B_freq_cgs( Temp_k, freq_value )
	B_freq_cgsLIST.append(B_freq_cgs)


#Normalized lists 
max_b_wave = max(B_wave_cgsLIST)
for idx3, B_value in enumerate(B_wave_cgsLIST):
	B_waveNORMLIST.append(B_value/max_b_wave)

max_b_freq = max(B_freq_cgsLIST)
for idx4, B_value in enumerate(B_freq_cgsLIST):
	B_freqNORMLIST.append(B_value/max_b_freq)


"""
=====================================================
Integrate Spectrum from Arbitrary Bounds
=====================================================
"""
#Sets the DEFAULT bounds for integration
print("Setting integration bounds to default... ")
low_bound_wave  = 1
high_bound_wave = 1
#Set to whatever you would like: gets executed by raising to power. 
print("Setting integration bounds to cover all values")
low_bound_wave  = wave_peak_log - 2
high_bound_wave = wave_peak_log + 2


"""
print("TEMP SET")
integration_lower_cm = calc_wavelength_given_energy(999.)
integration_upper_cm = calc_wavelength_given_energy(1001.)
low_bound_wave = np.log10(integration_lower_cm)
high_bound_wave = np.log10(integration_upper_cm)
"""

#low_bound_freq  = freq_peak_log - 1
#high_bound_freq = freq_peak_log + 1


integrated_intensity_cgs = integrate_basic( B_wave_cgsLIST, boundwaveLIST, low_bound_wave, high_bound_wave )
print("===================================================")
print("The Intensity (cgs) wave for [  %.1f , %.1f ] is : %g" % (low_bound_wave, high_bound_wave, integrated_intensity_cgs))
#integrated_intensity_freq_cgs = integrate_basic( B_freq_cgsLIST, boundfreqLIST, low_bound_freq, high_bound_freq )
#print("The Intensity (mks) freq for [  %.1f , %.1f ] is : %g" % (low_bound_freq, high_bound_freq, integrated_intensity_freq_cgs))
radius = 5. 
intensity_cgs_value =  calc_B_wavelength_cgs( Temp_k, specific_wave)
intensity_eyecgs_value =  calc_B_wavelength_cgs( Temp_k, specific_wave_eye)

#print("The Intensity (cgs) specific wave for [  %.3e ]cm is : %.5e" % ( specific_wave , intensity_cgs_value))
print("The Intensity (cgs) specific wave for [  %.3e ]cm is : %.5e" % ( specific_wave_eye , intensity_cgs_value))
print("The Intensity (mks) specific wave for [  %.3e ] m is : %.5e" % ( specific_wave/100. , calc_B_wavelength_cgs( Temp_k, specific_wave/100.)))
print("Specific Wavelength ==> %g Angstroms" % (specific_wave * (10**8.))) 


#Find Photons emitted from an eye. 
energy_photon_joules = h_cgs * c_cgs / specific_wave_eye 
energy_photon_ergs = energy_photon_joules * (10**7.)
eye_radius = 1. 
eye_surface_area = get_surfacearea_cgs( eye_radius )
num_eyephotons = intensity_eyecgs_value * eye_surface_area / energy_photon_ergs
print("Peak [cm]= %.7f " % specific_wave_eye)
print("The # of photons emitted from the eye at peak: %.4e " % num_eyephotons)
#Find Photons emitted from a potato. 
energy_photon_joules = h_cgs * c_cgs / specific_wave 
energy_photon_ergs = energy_photon_joules * (10**7.)
potato_radius = 5. 
potato_surface_area = get_surfacearea_cgs( potato_radius )
num_photons = intensity_cgs_value * potato_surface_area / energy_photon_ergs
print("Peak [cm]= %.7e " % specific_wave)
print("The # of photons emitted from a potato: %.4e " % num_photons)
#in cm: http://faculty.buffalostate.edu/sabatojs/courses/GES639/S10/reading/mass_luminosity.pdf
star_radius = 1.33*((25)**.555) * 69 * (10**9.)
star_surface_area = get_surfacearea_cgs( star_radius )
num_photons = intensity_cgs_value * star_surface_area / energy_photon_ergs
print("The # of photons emitted from a star: %.4e " % num_photons)


surfacearea_cgs = get_surfacearea_cgs( radius  )

#surfacearea_mks = get_surfacearea_mks( radius  )
print("The Intensity (cgs) given a SA [  %.1f ]cm^3 is : %g" % (surfacearea_cgs, integrated_intensity_cgs*surfacearea_cgs))
print("===================================================")





"""
----------------------------------------------
# BEGIN PLOT FUNCTION Calls. { X,  Y}
----------------------------------------------
"""
"""
#======================================================
# 1 A. Intensity per Wavelength MKS
#======================================================
title_label = "Specific Intensity vs. mks Wavelength: %s F " % str(Temperture_F)
x_label = "Wavelength [m]"
y_label = "Intensity [W sr-1 m-3]"
x_data  = boundwaveLIST
y_data  = B_wave_mksLIST
legend_val = 0
pointsize = 3
yflip = False
ylog = 0
xlog = 1

pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, \
         legend_val, pointsize, xlog, ylog,  yflip, pcount)
"""

#======================================================
# 2 A. Intensity per Wavelength CGS
#======================================================
title_label = "Specific Intensity vs. cgs Wavelength: %s F " % str(Temperture_F)
x_label = "Wavelength [cm]"
y_label = "Intensity [ergs s-1 sr-1 cm-3]"
x_data  = boundwaveLIST
y_data  = B_wave_cgsLIST
legend_val = 0
pointsize = 3
yflip = False
ylog = 0
xlog = 1

pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, \
         legend_val, pointsize, xlog, ylog,  yflip, pcount)



#======================================================
# A. Intensity per Wavelength Normalized
#======================================================
title_label = "Specific Intensity vs. Wavelength: %s F " % str(Temperture_F)
x_label = "Wavelength [cm]"
y_label = "Normalized Intensity [W sr-1 cm-3]"
x_data  = boundwaveLIST
y_data  = B_waveNORMLIST
legend_val = 0
pointsize = 3
yflip = False
ylog = 0
xlog = 1

pcount = plot_basic(x_data, y_data, title_label, x_label, y_label, \
         legend_val, pointsize, xlog, ylog,  yflip, pcount)

"""
#======================================================
# A. part II Intensity per Frequency MKS
#======================================================
title_label = "Specific Intensity vs. Freq: %s F"  % str(Temperture_F)

x_label = "Frequency [Hz]"
y_label = "Intensity [W sr-1 m-2 Hz-1]"
x_data  = boundfreqLIST
y_data  = B_freq_mksLIST
legend_val = 0
pointsize = 3
yflip = False
ylog = 0
xlog = 1

pcount = plot_basic(x_data, y_data, title_label, x_label, y_label,  \
	     legend_val, pointsize, xlog, ylog,  yflip, pcount)
"""

#======================================================
# A. part II Intensity per Frequency CGS
#======================================================
title_label = "Specific Intensity vs. Freq: %s F"  % str(Temperture_F)
x_label = "Frequency [Hz]"
y_label = "Intensity [erg s-1 sr-1 cm-2 Hz-1]"
x_data  = boundfreqLIST
y_data  = B_freq_cgsLIST
legend_val = 0
pointsize = 3
yflip = False
ylog = 0
xlog = 1

pcount = plot_basic(x_data, y_data, title_label, x_label, y_label,  \
	     legend_val, pointsize, xlog, ylog,  yflip, pcount)

#======================================================
# A. part II Intensity per Frequency Normalized
#======================================================
title_label = "Specific Intensity vs. Freq: %s F"  % str(Temperture_F)
x_label = "Frequency [Hz]"
y_label = "~ Normalized Intensity [W sr-1 m-2 Hz-1]"
x_data  = boundfreqLIST
y_data  = B_freqNORMLIST
legend_val = 0
pointsize = 3
yflip = False
ylog = 0
xlog = 1

pcount = plot_basic(x_data, y_data, title_label, x_label, y_label,  \
	     legend_val, pointsize, xlog, ylog,  yflip, pcount)





#-------------------------------------------
# Subsample -19
#-------------------------------------------
print("Time: ")
print (datetime.now() - startTime)
print("End of Program. ")






# Garbage comment for sublime convienence akjsdfjan


