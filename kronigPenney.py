#!/usr/bin/python

''' This script is used to analyze the energy levels predicted by the Kronig-
    Penney model of electron band structure for a periodic potential.  The end
	result of this program are plots of the E-K bands of the potential.
	
	/Kronig-Penney Model Potential/

	                          U
                              ^
                 b            |             a
              <---->          |          <----->
	----+     +----+     +----+     +----+     +----+     +---- Uo
        |     |    |     |    |     |    |     |    |     |
...     |     |    |     |    |     |    |     |    |     |    ...
        |     |    |     |    |     |    |     |    |     |
---------------------------------------------------------------> x
                      x = -b  |   x = a
                              |

	
	*** NOTE ***
	If 'Uo' is negative (instead of having ridges, there are pits), then the
	resulting energy banding is the same and can be obtained by swapping 'a'
	with 'b' and by shifting all energy values down by the magnitude of 'Uo'.

'''

from numpy import *
from matplotlib.pyplot import *

#------------------------------------------------------------------------------
# Physical Constants
#------------------------------------------------------------------------------
_c = 2.998e17			# Speed of light in nm/s
_m = 0.511e6/_c**2		# Electron mass in eV/c^2
_hbar = 6.583e-16		# h/(2*pi) in eV*s
_Uo = _hbar**2/(2*_m)	# If Uo = _Uo then u = 1

#------------------------------------------------------------------------------
# Model Parameters
# (units: a,b = nm, Uo,dE = eV)
#------------------------------------------------------------------------------
Uo = _Uo
a = pi
b = pi
Emax = 4*Uo
dE = 0.00001*Emax
gap_tol = 1.000001*dE

# Display model parameters
print "----------------------------------------------------------------------"
print "Model Parameters"
print "Uo (eV):", Uo
print "a (nm):", a
print "b (nm):", b
print "dE (eV):", dE
print "Energy gap tolerance (eV):", gap_tol

#------------------------------------------------------------------------------
# Function definitions
#------------------------------------------------------------------------------
def f_bound(E, Uo, a, b):
	u = sqrt(2*_m*Uo)/_hbar
	xi = E/Uo
	alpha = a*u*sqrt(xi)
	beta = b*u*sqrt(1 - xi)
	C = (1 - 2*xi)/(2*sqrt(xi*(1 - xi)))
	return cos(alpha)*cosh(beta) + C*sin(alpha)*sinh(beta)

def f_free(E, Uo, a, b):
	u = sqrt(2*_m*Uo)/_hbar
	xi = E/Uo
	alpha = a*u*sqrt(xi)
	beta = b*u*sqrt(xi - 1)
	C = (1 - 2*xi)/(2*sqrt(xi*(xi - 1)))
	return cos(alpha)*cos(beta) + C*sin(alpha)*sin(beta)

def f_top(Uo, a, b):
	u = sqrt(2*_m*Uo)/_hbar
	alpha = a*u
	C = -0.5*b*u
	return cos(alpha) + C*sin(alpha)

def f_bot(Uo, a, b):
	u = sqrt(2*_m*Uo)/_hbar
	beta = b*u
	C = 0.5*a*u
	return cosh(beta) + C*sinh(beta)

def mass_eff(k,E):
	return ((_hbar*_c)**2)/(polyfit(k,E,2)[0])


#==============================================================================
# Perform Calculations - The main part of the script
#==============================================================================

# Assign initial values of E
if Uo > 0:
	E = arange(0,Emax,dE)
elif Uo < 0:
	print "Uo must be greater than 0!"
	exit()
else:
	print "Uo cannont equal zero!"
	exit()


# Calculate the energy bands
f = []
E_bands = []
f_bands = []
for En in E:
	if En == 0:
		fn = f_bot(Uo, a, b)
		f.append(fn)
		if abs(fn) <= 1:
			E_bands.append(En)
			f_bands.append(fn)
	elif En < Uo:
		fn = f_bound(En, Uo, a, b)
		f.append(fn)
		if abs(fn) <= 1:
			E_bands.append(En)
			f_bands.append(fn)
	elif En == Uo:
		fn = f_top(Uo, a, b)
		f.append(fn)
		if abs(fn) <= 1:
			E_bands.append(En)
			f_bands.append(fn)
	elif En > Uo:
		fn = f_free(En, Uo, a, b)
		f.append(fn)
		if abs(fn) <= 1:
			E_bands.append(En)
			f_bands.append(fn)
f = array(f)				# f as a function of E
E = array(E)				# Emin to Emax by dE 
xi = E/Uo					# Emin/Uo to Emax/Uo by dE/Uo
E_bands = array(E_bands)	# Allowed Energies
f_bands = array(f_bands)	# Allowed values of f
xi_bands = E_bands/Uo		# Normalized allowed energies
K = arccos(f_bands)/(a+b)	# K values corresponding to energies
Kext = zeros(K.shape)		# Extended K values


# Find Band Gaps And Effective Mass At Band Ends

E_gaps = []					# Array of tuples containing band gaps
							# (Bottom Energy, Top Energy, Gap Width)

gap_idx = []				# Index of the bottom edge of every band

n = 0
for i in range(1,E_bands.size):
	chgE = E_bands[i] - E_bands[i-1]
	
	if chgE > gap_tol:
		E_gaps.append((E_bands[i-1], E_bands[i], chgE))
		gap_idx.append(i)
		n += 1 

	if n != 0:
		if n%2:
			Kext[i] = (n + 1)*pi/(a+b) - K[i]
		else:
			Kext[i] = n*pi/(a+b) + K[i]
	else:
		Kext[i] = K[i]
del n

# Normalized E_gaps
xi_gaps = [(gap[0]/Uo, gap[1]/Uo, gap[2]/Uo) for gap in E_gaps]

# Maximum value of K
Kmax = max(Kext)

# Energies of a free particle over the range -Kmax to Kmax
fakeK = arange(-Kmax,Kmax,0.001*Kmax)
E_empty = (_hbar*fakeK)**2/(2*_m)

# Calculate the effective masses at each band edge
effMass_gnd = mass_eff(Kext[0:3],E_bands[0:3]) # Ground State Effective Mass

effMass = []				# Array of tuples of effective masses for bands
							# [(mass at bottom, mass at top), (..),...]
							#  ----------"Band 1"-----------, ... 
mbot = effMass_gnd
mtop = 0
for i in gap_idx:
	mtop = mass_eff(Kext[i-3:i],E_bands[i-3:i])
	effMass.append((mbot,mtop))
	mbot = mass_eff(Kext[i:i+3],E_bands[i:i+3])
del mbot, mtop
#==============================================================================


#------------------------------------------------------------------------------
# Display Ground State Energy, Energy Gaps, & Effective masses at band ends
#------------------------------------------------------------------------------
print "----------------------------------------------------------------------"
print "Ground State Energy (E/Uo):", E_bands[0]/Uo
print "Energy Gaps (E/Uo)"
print "[Gap #]  (Lower Energy, Upper Energy, Energy Gap)"
for i in range(0,len(E_gaps)):
	print "[{0:2d}] ".format(i+1), xi_gaps[i]
print "----------------------------------------------------------------------"
print "Ground State Energy (eV):", E_bands[0]
print "Energy Gaps (eV)"
print "[Gap #]  (Lower Energy, Upper Energy, Energy Gap)"
for i in range(0,len(E_gaps)):
	print "[{0:2d}] ".format(i+1), E_gaps[i]
print "----------------------------------------------------------------------"
print "Ground State Effective mass (MeV):", effMass_gnd/1e6
print "Effective Mass At Band Edges (MeV)"
print "[Band #]  (Mass At Bottom, Mass At Top)"
for i in range(0,len(effMass)):
	print "[{0:2d}] ".format(i+1), (effMass[i][0]/1e6, effMass[i][1]/1e6)

#------------------------------------------------------------------------------
# Display Pretty Figures
#------------------------------------------------------------------------------
# Graphical Solution
'''
f1 = figure(1)
plot(xi, f, 'b',
	 xi, zeros(xi.shape), 'k',
	 xi, ones(xi.shape), 'r',
	 xi, -1*ones(xi.shape),'r',
	 figure=f1)
grid(True)
ylim([-4,4])
xlabel("E/Uo")
ylabel("f(E)")
'''

# Normalized E-K Plots
f2 = figure(2)
plot(-1*K*(a+b)/pi, xi_bands, 'b', K*(a+b)/pi, xi_bands, 'b',
	 -1*Kext*(a+b)/pi, xi_bands, 'r', Kext*(a+b)/pi, xi_bands, 'r',
	 fakeK*(a+b)/pi, E_empty/Uo, 'g--',
	 figure=f2)
grid(True)
xlim([-1.1*Kmax*(a+b)/pi,1.1*Kmax*(a+b)/pi])
ylim([0,Emax/Uo])
xlabel("K*(a + b)/Pi")
ylabel("E/Uo")
title("EK Plots")

'''
# E-K Plots
f3 = figure(3)
plot(-1*K, E_bands, 'b', K, E_bands, 'b',
	 -1*Kext, E_bands, 'r', Kext, E_bands, 'r',
	 fakeK, E_empty, 'g--',
	 figure=f3)
grid(True)
xlim([-1.1*Kmax, 1.1*Kmax])
ylim([0,Emax])
xlabel("K  (1/nm)")
ylabel("E  (eV)")
title("EK Plots")
'''

show()

