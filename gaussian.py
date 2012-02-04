#! /usr/bin/python
"""
Advanced example of a plot done using matplotlib

This example shows how to control many of the aspects of a figure/plot done with
the matplotlib module.  It can be used as a template for other plots.
"""
import sys, scipy as sp, matplotlib.pyplot as plt

# minimum  and maxium values
xMin = -5
xMax = 5
yMin = -1
yMax = 1
xRange = xMax - xMin
yRange = yMax - yMin

# change in x values
dx = 1e-5

# mean
u = 0

# standard deviation
s = 1

# Calculate x and y
#  y is a gaussian distribution
x = sp.arange(xMin, xMax, dx)
y = ((2*sp.pi*s**2)**-0.5)*sp.exp(-(x-u)**2/(2*s**2))

# Create a figure, set axis properties, turn on the grid, and hold the current
# state of the axes
xTickw = float(xRange)/10	# difference between x ticks
yTickw = float(yRange)/10	# difference between y ticks
fig_height = 10				# height in inches
fig_width = 12				# width in inches
fig = plt.figure(num=1,
				figsize=(fig_width, fig_height),
				dpi=80,
				facecolor='w',
				edgecolor='k')
ax = fig.gca()
ax.set_xticks(sp.arange(xMin, xMax+xTickw, xTickw))
ax.set_yticks(sp.arange(yMin, yMax+yTickw, yTickw))
# ax.set_xscale('log')
ax.grid(True)
ax.hold(True)

# Add the plot to the figure, set titles, and display it
plt.plot(x, y, "b", linewidth=2, figure=fig)
plt.axis([xMin, xMax, yMin, yMax])
plt.title('Probability Distribution\n'
		+ r'$\sigma = ' + str(s) + r'$ '
		+ r'$\mu = ' + str(u) + r'$',
		family='serif', size=20, weight='bold')
plt.xlabel('X Axis', family='serif', size=16, style='italic')
plt.ylabel('Y Axis', family='serif', size=16, style='italic')
plt.xticks(family='serif', size=14)
plt.yticks(family='serif', size=14)
plt.show()


