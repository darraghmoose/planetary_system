import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from orbits import initial_positions, initial_velocity, body_update

#constants
AU = 1.495978707e11
G = 6.67e-11
M_SUN = 1.98892e30
timestep = 12 * 3600 #3600 seconds = 1 hour
num_steps = 3 * (365 * 24) 

#constants & settings
body_names = np.array(["Sun", "Earth"])
masses = np.array([1.98892e30, 5.9742E24])
initial_radii = np.array([0, 1.495978707e11])
plot_colours = ['#ef0000', '#3c6690']

#create blank arrays of equal dimensions to store position and velocity data
x_positions = np.zeros((num_steps,len(body_names))); y_positions = np.zeros((num_steps,len(body_names)))
x_velocities = np.zeros((num_steps,len(body_names))); y_velocities = np.zeros((num_steps,len(body_names)))


#initialise position  and velocities of all bodies surrounding the sun
x_positions[0], y_positions[0] = initial_positions(initial_radii)
x_velocities[0], y_velocities[0] = initial_velocity(x_positions[0], y_positions[0], masses)

#initialise velocity of sun to give system zero net momentum
x_velocities[0, 0] = -np.sum(masses[1:] * x_velocities[0, 1:]) / masses[0]
y_velocities[0, 0] = -np.sum(masses[1:] * y_velocities[0, 1:]) / masses[0]

#main simulation loop
crossings = [[] for i in range(len(body_names)-1)]
for loop in tqdm(range(1, num_steps), ncols=90, ascii=" ░▒█", desc="Running Simulation", colour="#00bb82"):  #"░▒█" -- " ▖▘▝▗▚▞█"
    x_positions[loop], y_positions[loop], x_velocities[loop], y_velocities[loop] = body_update(x_positions[loop-1], y_positions[loop-1], x_velocities[loop-1], y_velocities[:][loop-1], timestep, masses)
    #check for axis crossing
    for planet in range(1, len(body_names)):
        y_coord_now = y_positions[loop][planet] - y_positions[loop][0]
        y_coor_prev = y_positions[loop-1][planet] - y_positions[loop-1][0]
        x_coord_now = x_positions[loop][planet] - x_positions[loop][0]

        #look for sign change across axis 
        if y_coor_prev * y_coord_now < 0 and x_coord_now > 0:
            #timestep at [loop] = loop * timesteps
            time_after = loop * timestep
            crossing_time = time_after - (y_coord_now)*(timestep / (y_coord_now - y_coor_prev))
            crossings[planet-1].append(crossing_time)


#print period data in table
print(f'Orbital period results:\n')
print(f'{'Planet':<12}{'Seconds':>15}{'Hours':>15}{'Days':>15}')
print(f'---------------------------------------------------------')
for i in range(1, len(body_names)):
    times = crossings[i-1]
    if len(times) < 2: #body needs to pass x axis twice to measure period
        print(f'{body_names[i]:<12}{'Insufficient crossings':>45}')
        continue
    orbital_period = np.diff(times).mean()
    print(f'{body_names[i]:<12}{orbital_period:>15,.1f}{orbital_period/3600:>15,.2f}{orbital_period/(3600*24):>15,.2f}')

#plot orbits of bodies to screen
plt.figure(figsize=(7, 7))
for i in range(len(body_names)):
    plt.plot(x_positions[:, i]/AU, y_positions[:, i]/AU, label=body_names[i], color=plot_colours[i])
plt.legend()
plt.axis('equal')
plt.xlabel('X (AU)', fontsize=16)
plt.ylabel('Y (AU)', fontsize=16)
"""
#create zoomed in plot of sun's orbit
axins = plt.inset_axes([0.6, 0.6, 0.37, 0.37])   #section of plot to zoom in on
axins.plot(x, y1, linestyle='-', label='Sample 1')
axins.plot(x, y2, linestyle='--', label='Sample 2')

sub_axes = plt.axes([0.3, 0.55, 0.15, 0.15]) #type: ignore
sub_axes.plot(x_positions[:,0], y_positions[:,0], color=plot_colours[0])
sub_axes.indicate_inset_zoom(sub_axes)
"""
plt.savefig("images/unstructured_1_moving_sun.pdf")
plt.show()