import math
import numpy as np
import matplotlib.pyplot as plt

#This will be the Sun and Earth only model
AU = 1.495978707e11 #meters
G = 6.67e-11 #N m^2 kg^-2
M_SUN = 1.98892e30 #kg
timestep = 1000 #seconds

#declare starting positions anf velocities of planetary bodies
x_sun = y_sun = 0
x_velocity_sun = y_velocity_sun = 0

x_earth = 1 * AU
y_earth = 0
x_velocity_earth = 0
y_velocity_earth = -1 * math.sqrt(G * M_SUN / AU ) #clockwise orbits

def acceleration(x_position, y_position):
    seperation = math.sqrt(x_position**2 + y_position**2)
    k = -G * M_SUN / seperation**3
    return x_position * k, y_position * k

def earth_update(x_position, y_position, x_velocity, y_velocity, timestep):
    """
    """
    #calculate all a values
    a_x_position , a_y_position = x_velocity,  y_velocity 
    a_x_velocity, a_y_velocity = acceleration(x_position, y_position)

    #calculate all b values
    b_x_position, b_y_position = x_velocity + 0.5*timestep*a_x_velocity, y_velocity + 0.5*timestep*a_y_velocity
    b_x_velocity, b_y_velocity = acceleration(x_position + 0.5*timestep*a_x_position, y_position + 0.5*timestep*a_y_position)

    #calculate all c values
    c_x_position, c_y_position = x_velocity + 0.5*timestep*b_x_velocity, y_velocity + 0.5*timestep*b_y_velocity
    c_x_velocity, c_y_velocity = acceleration(x_position + 0.5*timestep*b_x_position, y_position + 0.5*timestep*b_y_position)
   
    #calculate all d values
    d_x_position, d_y_position = x_velocity + timestep*c_x_velocity, y_velocity + timestep*c_y_velocity
    d_x_velocity, d_y_velocity = acceleration(x_position + timestep*c_x_position, y_position + timestep*c_y_position)
   

    #calculate weighted average of values to update all variables
    updated_x_position = x_position + ((timestep/6) * (a_x_position + 2*b_x_position + 2*c_x_position + d_x_position))
    updated_y_position = y_position + ((timestep/6) * (a_y_position + 2*b_y_position + 2*c_y_position + d_y_position))
    updated_x_velocity = x_velocity + ((timestep/6) * (a_x_velocity + 2*b_x_velocity + 2*c_x_velocity + d_x_velocity))
    updated_y_velocity = y_velocity + ((timestep/6) * (a_y_velocity + 2*b_y_velocity + 2*c_y_velocity + d_y_velocity))
    return updated_x_position, updated_y_position, updated_x_velocity, updated_y_velocity

x_plotting_earth = [x_earth/AU]
y_plotting_earth = [y_earth]
counter=0

for i in range(35000):
    x_earth, y_earth, x_velocity_earth, y_velocity_earth = earth_update(x_earth, y_earth, x_velocity_earth, y_velocity_earth, timestep)
    x_plotting_earth.append(x_earth/AU)
    y_plotting_earth.append(y_earth/AU)
    if -0.0001*AU<y_earth<0.0001*AU and 1-0.0001*AU<x_earth<1.0001*AU:
        print("success")
        year_time = i
        break

print(f'Year time: {i*timestep}')

print(counter)

plt.figure(figsize=(8, 8))
plt.plot(x_plotting_earth, y_plotting_earth)
plt.xlim(-1.05, 1.05)
plt.ylim(-1.05, 1.05)
plt.show()

"""
Want to find when 1 orbit has been complete










Problem: seperation is increasing by oders of magnitude and eventually blows up (E+104)
guess: issue is with update order (updating all 4 variables at once)
FIXED: doc outlining formulae has V * x_i when it should be v_{xi}. Fuming
good sign is position and velocity are remaining mostly constant across 20 iterations. 

problem: plotting x against y shows curve far from a circle
while magnitude of velocity seems to be accurate, direction seems wrong
FIXED: just needed more iterations
"""