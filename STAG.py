# Code to take dipole positions from a file created by CORAl and visualise them in 3D. 
#
# INSTRUCTIONS FOR USE:
# The input file "dipole_coords.txt" dipoles will be created automatically whenever CORAl analyses a DDSCAT-type 
# shape.dat file. Simply run this python script from the command line in the same folder to read this input file
# and visualise the shape in 3D. The plot is interactive and fully rotatable.
#
# This code also displays a Mie sphere of equivalent volume to the shape as a to-scale comparison.
#
# The code works on the mac it was designed on, but has not been tested on any other computers - please e-mail
# the following e-mail address if you have any questions or issues: matthew.g.lodge@gmail.com
#
# TROUBLESHOOTING
#
# This exact version is designed to produce two windows (side-by-side) on mac computers; however, this will not
# work if running on windows. Simply comment out any lines that produce errors as a workaround, e.g:
#
# - matplotlib.use("TkAgg") # use a backend to allow the current_fig_manager section to position the windows
# - plt.get_current_fig_manager().window.wm_geometry("+750+100")
#
# and then uncomment the lines under "windows version".
# 
# Matt Lodge 01/07/22

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib
import pandas as pd
matplotlib.use("TkAgg") # use a backend to allow the current_fig_manager section to position the windows
import matplotlib.pyplot as plt


print("\n\n ----------------------------------------------------------------------------------------")
print("            Welcome to S.T.A.G (Simulated Three-dimensional Aerosol Geometries!)")
print(" ----------------------------------------------------------------------------------------")


# Import dipole positions
dipoles = pd.read_csv('dipole_coords.txt', header=None, names=['X', 'Y', 'Z'])

print(" ",len(dipoles)-1," dipoles imported successfully.")


# The final row gives information about the lattice dimensions, the number of dipoles, and the distance between dipoles - store this info before moving on
STAG_lattice_dim=dipoles['X'][len(dipoles)-1] # store the lattice dimension value
N=dipoles['Y'][len(dipoles)-1] # store the number of dipoles N from our program
d=dipoles['Z'][len(dipoles)-1]/1000000 # store the value of the spacing between dipoles "d". This is usually a float, but we need to import it as an int as part of the array to keep it as one import file - we export it from our program truncated in pico-metre units (so we still have the data in micrometres to 6 dp) and when we import it here, divide it by 1000000 to turn back into um units. This truncates some info but is accurate enough to plot on a graph for visuals here.



# create a 3D grid composed entirely of zeroes, with our lattice dimensions
grid= np.zeros((STAG_lattice_dim,STAG_lattice_dim,STAG_lattice_dim),dtype=int) # create a grid composed of zeroes
print(" ",STAG_lattice_dim,"x",STAG_lattice_dim,"x",STAG_lattice_dim," grid created.\n")


# Now that we know the number of dipoles "N", scan through all rows of the imported integer matrix "dipoles" and change any points in our zero-array "grid" from 0->1 wherever there are dipole positions (at any dipole coordinate)
for i in range (N):
    grid[dipoles['X'][i]][dipoles['Y'][i]][dipoles['Z'][i]]=1  


# create a custom coordinate map to display coords and shift them to center around 0
cx,cy,cz = (np.indices((STAG_lattice_dim+1,STAG_lattice_dim+1,STAG_lattice_dim+1)) - (STAG_lattice_dim/2))*d   #the coords should have a max value of "lattice_dim+1", and then we want to re-center them "-(lattice_dim/2)" 

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

plt.get_current_fig_manager().window.wm_geometry("+750+100")

# set axis titles
ax.set_xlabel('X (\u03BCm)')
ax.set_ylabel('Y (\u03BCm)')
ax.set_zlabel('Z (\u03BCm)')

# plot voxels - the cx,cy,cz arguments are optional, to give the correct units on the axes. We could leave these out, and it would just give integer values from 0->lattice_dimension on each axis (e.g. 0-10)
ax.voxels(cx,cy,cz, grid, edgecolor="k")

ax.set_axis_off() # optional: removes axes and grey area around shape


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

plt.get_current_fig_manager().window.wm_geometry("+50+100")


# CODE TO DISPLAY A MIE SPHERE

a_eff=pow((3.0*N/(4.0*np.pi)),(1.0/3.0))*d; # calculate the effective radius if the dipoles were arranged in a sphere of the same total volume

print(" A Mie sphere of the same volume would have a radius of ",round(a_eff,3),"\u03BCm. Plotting graphs for both...")

# Make data - create x,y,z positions in spherical polars (with radius = a_eff)
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
a = a_eff * np.outer(np.cos(u), np.sin(v)) # this is just x,y,z in spherical polars (e.g. x= r*cos(u)*sin(v))
b = a_eff * np.outer(np.sin(u), np.sin(v))
c = a_eff * np.outer(np.ones(np.size(u)), np.cos(v))

# set axes labels for the sphere
ax.set_xlabel('X (\u03BCm)')
ax.set_ylabel('Y (\u03BCm)')
ax.set_zlabel('Z (\u03BCm)')

# Plot the surface
ax.plot_surface(a, b, c)

print("\n Images loaded.\n\n")

# Match the axes to the fractal plot so the size difference between them is clear

side_length= (STAG_lattice_dim/2.0)*d # find the maximum side length of the grid

plt.xlim(-side_length,side_length) # set this as the limits for the x and y axes
plt.ylim(-side_length,side_length)
ax.set_zlim(-side_length,side_length) # z axes needs slightly different syntax



ax.set_axis_off() # optional: removes axes and grey area around shape



plt.show() # display all graphs!

