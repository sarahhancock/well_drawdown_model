import numpy as np
import matplotlib.pyplot as plt
from scipy.special import exp1 
import matplotlib.animation as animation
plt.style.use('ggplot')

"""
PROGRAM DESCRIPTION

This program models the impact of transient groundwater flow
and pumping wells on the hydraulic head in a confined aquifer.

It starts with a 100 by 100 numpy array containing the initial
hydraulic head levels in the aquifer, and it takes user input 
for the aquifer properties, the number of wells, and the time
for the simulation to be run.

It uses the finite difference method to calculate the drawdown
due to transient water flow at each point in the grid for each 
timestep, and it uses the Theis solution (calculated using the scipy 
package)  to calculate the drawdown due to each pumping well at each
point in the grid for each timestep. Then, because of the principle 
of superposition, the drawdown values at each point are summed to get
the total change in hydraulic head at each timestep.

Finally, matplotlib's imshow function is used to create a pcolormesh
animation of the hydraulic heads every 20 time steps, and this 
animation is saved to the current directory.

Assumes valid user input (i.e. users input integer values when asked).
"""



"""
METHODS
"""


"""
Returns the drawdown due to transient flow after one time step
in the aquifer at each point in the grid as a numpy array with 
the same size as the grid.
"""
def get_transient_drawdown(S,T,dx,H,dt):
    
    H_new = np.copy(H)
    b = T*dt/(S*dx**2)
    a = 1-4*b
    
    H_new[1:99, 1:99] = a*H[1:99, 1:99] + \
                        b*(H[2:100, 1:99] + H[0:98, 1:99] + \
                           H[1:99, 2:100] + H[1:99, 0:98])
    
    transient_drawdown = H_new - H
    return transient_drawdown

"""
Returns the drawdown due to a particular well at a 
single point at a distance, r, from the well.
Used as a helper function for get_well_drawdowns.
"""
def get_well_drawdown(S, T, r, t, Q):
    u = S*r**2/(4*T*(t+1**-9))
    W = exp1(u)
    return Q*W/(4*np.pi*T)


"""
Adds a well to the grid.
"""
def add_well(well_locations, i, j, Q, X):
    well_drawdown_prev = X*0
    well_locations.append([i,j,Q,well_drawdown_prev])
    return

"""
Gets the distance of a point from a well location.
Used as a helper function for get_well_drawdowns.
"""
def get_r_from_well(i,j,i_well, j_well):
    return np.sqrt((i-i_well)**2 + (j - j_well)**2)


"""
Finds the drawdown caused by one well at each point 
on the grid and returns the drawdown values as a numpy
array with the same shape as the grid.
"""
def get_well_drawdowns(well, S, T, t, dx, dy):
    i_well = well[0]
    j_well = well[1]
    Q = well[2]
    well_drawdown_prev = well[3]
    i = np.arange(0,100,1)
    j = np.arange(0,100,1)
    I, J = np.meshgrid(i,j)
    radii = np.copy(I)
    radii[:] = get_r_from_well(I*dx,J*dy,i_well, j_well)
    radii[int(j_well/dx), int(i_well/dy)] = 1**(-10)
    well_drawdown = np.copy(radii)
    well_drawdown = get_well_drawdown(S,T,radii[:],t,Q)
    added_well_drawdown = well_drawdown - well_drawdown_prev
    well_drawdown_prev[:] = well_drawdown
    well[3] = well_drawdown_prev
    return added_well_drawdown



if __name__ == '__main__':

    """
    This part of the code sets up the initial conditions for the hydraulic head in a 
    100 by 100 numpy grid. The dimensions of the aquifer are set here as 1000 by 1000,
    and the head has been raised between x = (200,300) and y = (600,800) to demonstrate
    the transient groundwater flow, but these initial conditions can easily be changed
    for different scenarios,
    """
    dx = 10
    dy = 10
    xLen = 1000
    yLen = 1000
    x = np.arange(0., xLen, dx)
    y = np.arange(0., yLen, dy)
    X,Y = np.meshgrid(x,y)
    
    H = np.full((100,100), 50.)
    H[20:30, 60:80] = 60.
    
    well_locations = []
    min_h = np.amin(H)
    max_h = np.amax(H)
    t = 0
    count = 0
    heads_over_time = []
    

    """
    This part of the code prompts the user to enter aquifer properties, well locations, 
    well properties, and the time for which the simulation should be run. It is important 
    to note that storativity values are typically very small (0.0001 to 0.00001), but this
    creates a very small time step, dt, for using the finite difference method to get the drawdown
    for transient flow. For testing the program, I recommend using a storativity of around 0.1 or
    a very small time over which the simulation is run. However, if this program was being used
    to model an actual aquifer over time with small storativity values, it would still work, it
    would just take a long time to run.
    """
    
    print("\nThis program allows you to simulate the effect of wells on a confined aquifer over time.  \n ")
    print("Your aquifer is {} by {} m. \n ".format(xLen,yLen))
    print("Now, you can input the properties of the aquifer.   \n")
    S = float(input("Enter the storativity (unitless):  "))
    T = float(input("Enter the transmissivity (m^2/day):  "))
      
    dt = S*dx**2/(4*T)
    
    n_wells = input("How many wells would you like to add?  ") 
    for i in range(0, int(n_wells)):
        x_well = int(input("Please enter the x coordinate of well " + str(i+1) + " in meters.  "))
        y_well = int(input("Please enter the y coordinate of well " + str(i+1) +  " in meters.  "))
        Q = float(input("Please enter Q, the pumping rate of the well, in m^3/day (Q < 0 for pumping water out, and Q > 0 for pumping water in).  "))
        add_well(well_locations, x_well, y_well, Q, X)
        
    t_max = float(input("Enter the number of days you would like the simulation to run for: "))
    print("Ok, running simulation for " + str(t_max) + " days.")
    
    
    while (t <= t_max+1):
        transient_drawdown = get_transient_drawdown(S,T,dx,H,dt)    
        for well in well_locations:
            well_drawdown = get_well_drawdowns(well, S, T, t, dx, dy)
            H += well_drawdown
        H += transient_drawdown
        H = H*(H>=0)
        
        if (count%20 == 0):
            heads_over_time.append([t, H])
            if np.amin(H) < min_h:
                min_h = np.amin(H)
            if np.amax(H) > max_h:
                max_h = np.amax(H)
                
        t+=dt
        count+=1
        
"""
This part of the code creates a pcolormesh animation of the hydraulic 
head over the given time range by creating a 3D array indexed by the 
x-coordinate, y-coordinate, and time.
"""

nPlots = len(heads_over_time)
H3D = np.full((100,100,nPlots),0.)
T = np.arange(0, t_max, dt)


for k in range(nPlots):
    heads = heads_over_time[k][1]
    for i in range(0,100):
        for j in range(0,100):
            H3D[i,j,k] = heads[i,j]

fig = plt.figure(facecolor='w', figsize = (10,6))
ax = plt.axes(xlabel = "X Position (m)", ylabel = "Y Position (m)")
ax.set_title("Hydraulic Head from t=0 to t={} days".format(int(t_max)))
im = ax.imshow(H3D[:,:,0], origin = 'lower', vmin = min_h, vmax = max_h, extent = (0, xLen, 0, yLen))
im.set_data(H3D[:,:,0])
plt.close()
fig.colorbar(im, label = "Hydraulic Head (m)")

def update_data(i):
    im.set_data(H3D[:,:,i])
    return

ani = animation.FuncAnimation(fig, update_data, interval = 150, frames = nPlots, blit = False)
filename = 'aquifer_model.mp4'
ani.save(filename)
print("Task completed. File saved as {} in current directory".format(filename))



"""
Sarah Hancock, seh2209
12/4/2019

***ADDITIONAL SOURCES (other than class materials)***

Calculating Well Drawdown Resources:

Introduction to Hydrology, Margulis (textbook) (pg. 303 especially)
    (for Theis solution equation and basic hydrology resources)
    
https://scipython.com/blog/linear-and-non-linear-fitting-of-the-theis-equation/)
    (for numerical solution of Theis solution in Python)

Drawdown Superposition Resources:

https://pubs.usgs.gov/of/1984/0459/report.pdf 

http://inside.mines.edu/~epoeter/_GW/15wh4Superpsition/WellHydraulics4pdf.pdf
    
Animating pcolormesh Resources:

https://brushingupscience.com/2016/06/21/matplotlib-animations-the-easy-way/

https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.imshow.html

"""