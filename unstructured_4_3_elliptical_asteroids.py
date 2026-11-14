import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from orbits import initial_positions, initial_velocity, body_update, elliptical_initial_velocity

def kinetic_energy(body_mass, body_x_velocity, body_y_velocity):
    energy = 0.5 * body_mass * (body_x_velocity**2 + body_y_velocity**2)
    return energy

#constants
AU = 1.495978707e11
G = 6.67e-11
M_SUN = 1.98892e30

#simulation runtime parameters
timestep = 24 * 3600        #seconds per timestep
num_steps = 300 * (365)  #iterations per simulaiton

#use arrays to store data for each body      -      data taken from https://solarsystem.nasa.gov/planet-compare/ 
body_names = np.array(["Sun", "Jupiter", "Asteroid"])
masses = np.array([1.98892e30, 1.89813E27, 0.95E21])
Jupiter_radius = 7.78479E11
kirkwood_radii = [2.5*AU, 2.83*AU, 2.96*AU, 3.28*AU]
eccentricities = np.array([1, 1, 0.95])

def main_simulation(asteroid_radius: float):
    """Simulates the Sun, Jupiter and a single asteroid starting at asteroid_radius.
    Returns the x and y positions of all bodies and the kinetic energy of the asteroid."""

    initial_radii = np.array([0, Jupiter_radius, asteroid_radius])

    #create blank arrays of equal dimensions to store position and velocity data
    x_positions = np.zeros((num_steps,len(body_names))); y_positions = np.zeros((num_steps,len(body_names)))
    x_velocities = np.zeros((num_steps,len(body_names))); y_velocities = np.zeros((num_steps,len(body_names)))

    #initialise position and velocity of all bodies surrounding the sun
    x_positions[0], y_positions[0] = initial_positions(initial_radii)
    x_velocities[0], y_velocities[0] = elliptical_initial_velocity(x_positions[0], y_positions[0], eccentricities)

    #initialise velocity of sun to give system zero net momentum
    x_velocities[0, 0] = -np.sum(masses[1:] * x_velocities[0, 1:]) / masses[0]
    y_velocities[0, 0] = -np.sum(masses[1:] * y_velocities[0, 1:]) / masses[0]

    #track KE of asteroid
    asteroid_energy = np.zeros(num_steps)
    asteroid_energy[0] = kinetic_energy(masses[2], x_velocities[0][2], y_velocities[0][2])

    #main simulation loop
    for loop in tqdm(range(1, num_steps), ncols=90, ascii=" ░▒█", desc=f"Asteroid at {asteroid_radius/AU:.2f} AU", colour="#00bb82"):  #" ░▒█" -- " ▖▘▝▗▚▞█"
        x_positions[loop], y_positions[loop], x_velocities[loop], y_velocities[loop] = body_update(x_positions[loop-1], y_positions[loop-1], x_velocities[loop-1], y_velocities[loop-1], timestep, masses)
        #Measure astroid KE to gague orbit stability
        asteroid_energy[loop] = kinetic_energy(masses[2], x_velocities[loop][2], y_velocities[loop][2])

    return x_positions, y_positions, asteroid_energy


#run simulation
all_x_positions = []; all_y_positions = []; all_asteroid_energy = []
for radius in kirkwood_radii:
    x_positions, y_positions, asteroid_energy = main_simulation(radius)
    all_x_positions.append(x_positions); all_y_positions.append(y_positions)
    all_asteroid_energy.append(asteroid_energy)

#print energy deviation data in table
normalised_energy = [energy/np.mean(energy) for energy in all_asteroid_energy]
deviation = [np.max(energy)-1 for energy in normalised_energy]

print(f'\nEnergy Deviation Data:\n')
print(f'{'Radius (AU)':<15}{'Deviation (%)':>15}')
print(f'------------------------------')
for i in range(len(kirkwood_radii)):
    print(f'{kirkwood_radii[i]/AU:<15,.2f}{deviation[i]*100:<15,.1f}')


#plot orbits of bodies to screen
plot_colours = ['#ef0000', '#c98226', '#8c8987', '#dfb18c', '#3c6690', '#ff6b4d', '#f8d280', '#5d8b8e', '#4e5a74']
plt.figure(1, figsize=(6, 6))
#sun and jupiter are near identical every run, so plot them from the first
plt.plot(all_x_positions[0][:, 0]/AU, all_y_positions[0][:, 0]/AU, color=plot_colours[0], label=body_names[0])
plt.plot(all_x_positions[0][:, 1]/AU, all_y_positions[0][:, 1]/AU, color=plot_colours[1], label=body_names[1])
for i in range(len(kirkwood_radii)):
    plt.plot(all_x_positions[i][:, 2], all_y_positions[i][:, 2], color=plot_colours[i+2], label=f'Asteroid {kirkwood_radii[i]/AU:.2f} AU')
plt.xlabel('X [AU]')
plt.ylabel('Y [AU]')
plt.legend()
plt.tight_layout()

#create normalised plot of asteroid energy
plt.figure(2, figsize=(6, 6))
plt.title("Asteroid Kinetic Energy")
for i in range(len(kirkwood_radii)):
    plt.plot(np.linspace(0, timestep*num_steps, num=num_steps), normalised_energy[i], color=plot_colours[i+2], label=f'Asteroid {kirkwood_radii[i]/AU:.2f} AU') 
plt.ylabel('Normalised Kinetic Energy')
plt.xlabel('Time')
plt.legend()
plt.tight_layout()

plt.show()