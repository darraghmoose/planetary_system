import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from orbits import initial_positions, initial_velocity, body_update

#constants & parameters
AU = 1.495978707e11
G = 6.67e-11
M_SUN = 1.98892e30
timestep = 12 * 3600 #3600 seconds = 1 hour
num_steps = 5 * (365 * 2)  #multiple of earth years

#Orbiting bodies data
body_names = np.array(["Sun", "Earth"])
masses = np.array([1.98892e30, 5.9742E24])
initial_radii = np.array([0, 1.495978707e11])
plot_colours = ['#ef0000', '#3c6690']

#create blank arrays of equal dimensions to store position and velocity data
x_positions = np.zeros((num_steps,len(body_names))); y_positions = np.zeros((num_steps,len(body_names)))
x_velocities = np.zeros((num_steps,len(body_names))); y_velocities = np.zeros((num_steps,len(body_names)))

#initialise position and velocity of all bodies surrounding the sun
x_positions[0], y_positions[0] = initial_positions(initial_radii)

#initialise velocities of all planets 
x_velocities[0], y_velocities[0] = initial_velocity(x_positions[0], y_positions[0], masses)


#main simulation loop
crossings = [[] for i in range(len(body_names)-1)]
for loop in tqdm(range(1, num_steps), ncols=90, ascii=" ░▒█", desc="Running Simulation", colour="#00bb82"):  #"░▒█" -- " ▖▘▝▗▚▞█"
    x_positions[loop], y_positions[loop], x_velocities[loop], y_velocities[loop] = body_update(x_positions[loop-1], y_positions[loop-1], x_velocities[loop-1], y_velocities[:][loop-1], timestep, masses)


#plot orbits of bodies to screen
plt.figure(figsize=(6, 6))
for i in range(len(body_names)):
    if i == 0: #need to use marker for sun as it never moves
        plt.plot(x_positions[0, i]/AU, y_positions[0, i]/AU, 'o', color=plot_colours[i], markersize=8, label=body_names[i]) #[0, i] used to avoid drawing same point over and over
    else:
        plt.plot(x_positions[:, i]/AU, y_positions[:, i]/AU, color=plot_colours[i], label=body_names[i])
plt.axis('equal')
plt.xlabel('X (AU)')
plt.ylabel('Y (AU)')
plt.legend()
plt.savefig("images/structured1.pdf")
plt.show()