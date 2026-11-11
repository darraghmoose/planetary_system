import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from orbits import initial_positions, initial_velocity, body_update

def kinetic_energy(body_mass, body_x_velocity, body_y_velocity):
    energy = 0.5 * body_mass * (body_x_velocity**2 + body_y_velocity**2)
    return energy

#constants
AU = 1.495978707e11
G = 6.67e-11
M_SUN = 1.98892e30

#simulation runtime parameters
timestep = 24 * 3600        #seconds per timestep
num_steps = 5 * (365 * 24)  #iterations per simulaiton

#use arrays to store data for each body      -      data taken from https://solarsystem.nasa.gov/planet-compare/ 
"""
body_names = np.array(["Sun", "Jupiter", "5:1 Resonance", "4:1 Resonance", "3:1 Resonance", "Asteroid Belt"])
masses = np.array([1.98892e30, 1.89813E27, 0.95E21, 0.95E21, 95E21, 95E21])
initial_radii = np.array([0, 5.204*AU, 1.780*AU, 2.065*AU, 2.502*AU, 2.8*AU])

body_names = np.array(["Sun", "Jupiter", "Unstable 1", "Unstable 2", "Asteroid Belt"])
masses = np.array([1.98892e30, 1.89813E27, 0.95E21, 0.95E21, 95E21])
initial_radii = np.array([0, 5.2*AU, (3.28 * AU), (4 * AU), 2.8*AU])       #initial_radii = np.array([0, 5.2*AU, (3.276 * AU), (4 * AU), 2.8*AU])

#Gravity assist
body_names = np.array(["Sun", "Jupiter", "Asteroid"])
masses = np.array([1.98892e30, 1.89813E27, 0.95E21])
initial_radii = np.array([0, 5.2*AU, 5.00* AU])  
"""

body_names = np.array(["Sun", "Jupiter", "Asteroid"])
masses = np.array([1.98892e30, 1.89813E27, 0.95E21])
initial_radii = np.array([0, 5.2*AU, 4.8*AU])


#create blank arrays of equal dimensions to store position and velocity data
x_positions = np.zeros((num_steps,len(body_names))); y_positions = np.zeros((num_steps,len(body_names)))
x_velocities = np.zeros((num_steps,len(body_names))); y_velocities = np.zeros((num_steps,len(body_names)))

#initialise position and velocity of all bodies surrounding the sun
x_positions[0], y_positions[0] = initial_positions(initial_radii)
x_velocities[0], y_velocities[0] = initial_velocity(x_positions[0], y_positions[0], masses)

#initialise velocity of sun to give system zero net momentum
x_velocities[0, 0] = -np.sum(masses[1:] * x_velocities[0, 1:]) / masses[0]
y_velocities[0, 0] = -np.sum(masses[1:] * y_velocities[0, 1:]) / masses[0]

#track KE of asteroids
asteroid_energy = np.zeros((num_steps, len(body_names)-2))
asteroid_energy[0] = kinetic_energy(masses[2:], x_velocities[0][2:], y_velocities[0][2:])

#main simulation loop
crossings = [[] for i in range(len(body_names)-1)]
for loop in tqdm(range(1, num_steps), ncols=90, ascii=" ░▒█", desc="Asteroid Simulation", colour="#00bb82"):  #" ░▒█" -- " ▖▘▝▗▚▞█"
    x_positions[loop], y_positions[loop], x_velocities[loop], y_velocities[loop] = body_update(x_positions[loop-1], y_positions[loop-1], x_velocities[loop-1], y_velocities[:][loop-1], timestep, masses)
    #Measure astroid KE to gague orbit stability
    asteroid_energy[loop] = kinetic_energy(masses[2:len(masses)], x_velocities[loop][2:len(masses)], y_velocities[loop][2:len(masses)])

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
normalised_energy = asteroid_energy/asteroid_energy[0]
deviation = np.max(normalised_energy, axis=0)

print(f'Orbital period results:\n')
print(f'{'Planet':<15}{'Seconds':>15}{'Hours':>15}{'Days':>15}{'Resonance':>15}{'Deviation (%)':>15}')
print(f'------------------------------------------------------------------------------------')
for i in range(1, len(body_names)):
    times = crossings[i-1]
    if len(times) < 2: #body needs to pass x axis twice to measure period
        print(f'{body_names[i]:<18}{'Insufficient crossings':>75}')
        continue
    orbital_period = np.diff(times).mean()
    if i == 1:
        jupiter_period = orbital_period
        print(f'{body_names[i]:<15}{orbital_period:>15,.1f}{orbital_period/3600:>15,.2f}{orbital_period/(3600*24):>15,.2f}{jupiter_period/orbital_period:>15,.3f}{'NA':>15}')
    else:
        print(f'{body_names[i]:<15}{orbital_period:>15,.1f}{orbital_period/3600:>15,.2f}{orbital_period/(3600*24):>15,.2f}{jupiter_period/orbital_period:>15,.3f}{abs(1-deviation[i-2])*100:15,.1f}')


#plot orbits of bodies to screen
plot_colours = ['#ef0000', '#c98226', '#8c8987', '#dfb18c', '#3c6690', '#ff6b4d', '#f8d280', '#5d8b8e', '#4e5a74']
plt.figure(1, figsize=(6, 6))
for i in range(len(body_names)):
    plt.plot(x_positions[:, i], y_positions[:, i], color=plot_colours[i], label=body_names[i])
plt.xlabel('X [m]')
plt.ylabel('Y [m]')
plt.legend()

#create normalised plot of asteroid energy
plt.figure(2, figsize=(6, 6))
plt.title("Asteroid Kinetic Energy")
for asteroid in range(len(normalised_energy[0])):
    plt.plot(np.linspace(0, timestep*num_steps, num=num_steps), normalised_energy[:, asteroid], color=plot_colours[asteroid+2], label=body_names[asteroid+2]) 
plt.ylabel('Normalised Kinetic Energy')
plt.xlabel('Time')
plt.legend()

#display plots
plt.show()