#! /usr/bin/python
'''
Creates interference patterns and displays them in an image

This program simulates point-sources that emit spherical radiation. Each source
in emits a spherical wave oscillating perpendicular to the source array plane.
The waves from each source are added up at each point in the image.  The
intensity is then calculated as the square of the summed radiation at each
point.  An image is then produced based upon the calculated intensities.
'''

import scipy as sp				# Needed for array objects
import numpy as np
import matplotlib.pyplot as plt	# Used to create the image


# Wavelength of radiation
wl = 1

# The real height and width of the image
height = 1e1
width = 1e1

# The configuration position vectors [x,y]
# (0,0 = top left corner of image)
d = 1.5

#n = int(height/d)
#pos = np.zeros([n,2])
#for i in range(1,n):
#	pos[i,0] = 0
#	pos[i,1] = d*i

ctr = (width/2, height/2)
pos = np.array([[ctr[0], ctr[1]],
				[ctr[0] + d, ctr[1] + d],
				[ctr[0] - d, ctr[1] + d],
				[ctr[0] + d, ctr[1] - d],
				[ctr[0] - d, ctr[1] - d]])



# Calculate important variables for computations
k = 2*sp.pi/wl
R = (height**2 + width**2)**0.5

# The size of a pixel
res = R*0.002

# Creates Amplitude Diffraction Pattern
img = np.zeros((round(height/res), round(width/res)))
for p in pos:
	for y in range(0,img.shape[0]):
		for x in range(0,img.shape[1]):
			# Calculates r
			r = ((x*res - p[0])**2 + (y*res - p[1])**2)**0.5

			# Calculate damping coefficient
			b = 1 - (r/R)**0.5

			# Calculates amplitude of the radiation field at current position
			img[y,x] += sp.cos(k*r)#*b


# Calculates Intensity Diffraction Pattern
img = img**2

# Normalize the intensity pattern
img /= np.max(img)

# Calculate image statistics
img_stats = {'max': np.max(img),
			 'min': np.min(img),
			 'mean': np.mean(img),
			 'std': np.std(img)}

# Print image statistics
print \
"""
min:		{0[min]}
max:		{0[max]}
mean:		{0[mean]}
std:		{0[std]}
""".format(img_stats)

# Displays Image
plt.imshow(img, cmap=plt.cm.gray)
plt.show()

