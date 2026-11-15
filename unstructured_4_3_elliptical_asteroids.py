import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from orbits import initial_positions, body_update, elliptical_initial_velocity, orbital_elements

#constants
AU = 1.495978707e11
G = 6.67e-11
M_SUN = 1.98892e30
YEAR = 365.25 * 24 * 3600

#simulation runtime parameters
timestep = 10 * 24 * 3600                   #seconds per timestep
num_steps =  int(20000 * 365.25 / 10)       #iterations per simulaiton
save_every = 500 

#use arrays to store data for each body      -      data taken from https://solarsystem.nasa.gov/planet-compare/ 
masses = np.array([1.98892e30, 1.89813E27, 0.95E21])
Jupiter_radius = 7.78479E11
kirkwood_radii = [2.5*AU, 2.83*AU, 2.96*AU, 3.28*AU, 2.6*AU] #2.6 is expected to be stable
eccentricities = np.array([0, 0.048, 0.15])

def main_simulation(asteroid_radius: float):
    """Simulates the Sun, Jupiter and a single asteroid starting at asteroid_radius.
    Returns the time, sun and asteroid positions, and the asteroid's semi-major axis and eccentricity."""

    np.random.seed(40451054) #helpful for debugging as get same random simulation each time

    initial_radii = np.array([0, Jupiter_radius, asteroid_radius]) * (1 + eccentricities)

    #current position and velocity of all bodies (only a few saved steps are stored)
    x_positions, y_positions = initial_positions(initial_radii)
    x_velocities, y_velocities = elliptical_initial_velocity(x_positions, y_positions, eccentricities)

    #initialise velocity of sun to give system zero net momentum
    x_velocities[0] = -np.sum(masses[1:] * x_velocities[1:]) / masses[0]
    y_velocities[0] = -np.sum(masses[1:] * y_velocities[1:]) / masses[0]
    
    #lists to store saved data
    times = []; semimajor_axis = []; eccentricity = []
    
    #main simulation loop
    for loop in tqdm(range(num_steps), ncols=90, ascii=" ░▒█", desc=f"Asteroid at {asteroid_radius/AU:.2f} AU", colour="#00bb82"):  #" ░▒█" -- " ▖▘▝▗▚▞█"
        if loop % save_every == 0:
            #Measure asteroid orbital elements 
            a, e = orbital_elements(x_positions[2], y_positions[2], x_velocities[2], y_velocities[2], x_positions[0], y_positions[0], x_velocities[0], y_velocities[0])
            times.append(loop * timestep / YEAR)
            semimajor_axis.append(a/AU); eccentricity.append(e)
        x_positions, y_positions, x_velocities, y_velocities = body_update(x_positions, y_positions, x_velocities, y_velocities, timestep, masses)
 
    return np.array(times), np.array(semimajor_axis), np.array(eccentricity)

#run simulation
all_times = []; all_a = []; all_e = []
for radius in kirkwood_radii:
    t, a, e = main_simulation(radius)
    all_a.append(a)
    all_e.append(e)
    all_times.append(t)

#print eccentricity data in table
print(f'\nEccentricity Data:\n')
print(f"{'Radius (AU)':<15}{'Max e':>10}{'Closest approach to Sun (AU)':>15}")
print(f'----------------------------------------')
for i in range(len(kirkwood_radii)):
    perihelion = all_a[i] * (1 - all_e[i])
    print(f'{kirkwood_radii[i]/AU:<15,.2f}{np.max(all_e[i]):>10.3f}{np.min(perihelion):>15.2f}') #if perihelion <1.7AU asteroid would have been collected by mars


#create plot of asteroid eccentricity against time
plot_colours = ['#dfb18c', '#ef0000', '#c98226',  '#4e5a74', '#8c8987', '#3c6690', '#ff6b4d', '#f8d280', '#5d8b8e']
plot_labels = ['3:1', '5:2', '7:3', '2:1', 'Stable']
plt.figure(figsize=(6, 6))
plt.title("Asteroid Eccentricity")
for i in range(len(kirkwood_radii)):
    plt.plot(all_times[i], all_e[i], color=plot_colours[i+2], label=f'{plot_labels}') 
plt.ylabel('Eccentricity', fontsize=16)
plt.xlabel('Time [years]', fontsize=16)
plt.legend(loc="upper left")
plt.tight_layout()
plt.savefig('images/unstructured_4_3_elliptical_eccentricities.pdf')
 
plt.show()