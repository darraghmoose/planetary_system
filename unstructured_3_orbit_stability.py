import numpy as np
import math
import matplotlib.pyplot as plt
from tqdm import tqdm
from orbits import initial_positions, initial_velocity, body_update

def kinetic_energy(body_mass: np.ndarray, body_x_velocity: np.ndarray, body_y_velocity: np.ndarray) -> np.ndarray:
    energy = 0.5 * body_mass * (body_x_velocity**2 + body_y_velocity**2)
    return energy 

#constants & settings
AU = 1.495978707e11
G = 6.67e-11
body_names = np.array(["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"])
masses = np.array([1.98892e30, 3.30104E23 , 4.86732E24, 5.9742E24, 6.4169E23, 1.89813E27, 5.68319E26, 8.68103E25, 1.0241E26])
initial_radii = np.array([0, 5.7909227E10 , 1.08209475E11, 1.495978707e11, 2.27956E11, 7.78340821E11, 1.426666422E12, 2.870658186E12, 4.498396441E12])
plot_colours = ['#ef0000', '#8c8987', '#dfb18c', '#3c6690', '#ff6b4d', '#c98226', '#f8d280', '#5d8b8e', '#4e5a74', '#5d8b8e']

timestep = 12 * 3600 #simulation time in seconds per timestep
num_steps = 4 * (365 * 24)

#create blank arrays of equal dimensions to store position and velocity data
x_positions = np.zeros((num_steps,len(body_names))); y_positions = np.zeros((num_steps,len(body_names)))
x_velocities = np.zeros((num_steps,len(body_names))); y_velocities = np.zeros((num_steps,len(body_names)))

#initialise position  and velocities of all bodies surrounding the sun
x_positions[0], y_positions[0] = initial_positions(initial_radii)
x_velocities[0], y_velocities[0] = initial_velocity(x_positions[0], y_positions[0], masses)

#initialise velocity of sun to give system zero net momentum
x_velocities[0, 0] = -np.sum(masses[1:] * x_velocities[0, 1:]) / masses[0]
y_velocities[0, 0] = -np.sum(masses[1:] * y_velocities[0, 1:]) / masses[0]

#track Kinetic energy of all planets around the Sun
planet_kinetic_energy = np.zeros((num_steps, len(body_names)-1))
#populate arrays with initial KE of all planets
planet_kinetic_energy[0] = kinetic_energy(masses[1:], x_velocities[0][1:], y_velocities[0][1:])

#main simulation loop
crossings = [[] for i in range(len(body_names)-1)]
for loop in tqdm(range(1, num_steps), ncols=90, ascii=" ░▒█", desc="Running Simulation", colour="#00bb82"):  #"░▒█" -- " ▖▘▝▗▚▞█"
    #update position and velocity of all bodies
    x_positions[loop], y_positions[loop], x_velocities[loop], y_velocities[loop] = body_update(x_positions[loop-1], y_positions[loop-1], x_velocities[loop-1], y_velocities[loop-1], timestep, masses)
    
    #calculate KE of all bodies
    planet_kinetic_energy[loop] = kinetic_energy(masses[1:len(masses)], x_velocities[loop][1:len(masses)], y_velocities[loop][1:len(masses)])
    
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
print(f'{'Planet':<12}{'Seconds':>15}{'Hours':>15}{'Days':>15}{'KE Deviation(%)':>18}')
print(f'---------------------------------------------------------------------------')
mean_planet_kinetic_energy = np.mean(planet_kinetic_energy, axis=0) #mean of column (planet) in array
normalised_planet_energy = planet_kinetic_energy/mean_planet_kinetic_energy
for i in range(1, len(body_names)):
    times = crossings[i-1]
    if len(times) < 2: #body needs to pass x axis twice to measure period
        print(f'{body_names[i]:<12}{'Insufficient crossings':>45}')
        continue
    orbital_period = np.diff(times).mean()
    
    #print(normalised_planet_energy[i])
    greatest_deviation = np.max(abs(normalised_planet_energy-1), axis=0)
    print(f'{body_names[i]:<12}{orbital_period:>15,.1f}{orbital_period/3600:>15,.2f}{orbital_period/(3600*24):>15,.2f}{greatest_deviation[i-1]*100:>18,.2f}')


#for very large simulations check if plotting is wanted
show_plot = 'y'
if num_steps > 5000:
    show_plot = input('Create plot? [y/n]: ') 

if show_plot == 'y':
    # #plot orbits of bodies
    plt.figure(1, figsize=(7, 7))
    for i in range(len(body_names)):
        plt.plot(x_positions[:, i]/AU, y_positions[:, i]/AU, label=body_names[i], color=plot_colours[i])
    plt.legend()
    plt.axis('equal')
    plt.xlabel('X (AU)')
    plt.ylabel('Y (AU)')
    plt.savefig("images/unstructured_3_total_system_studying_stability.pdf")

    #create normalised plot of planet energy
    #terrestrial planets
    plt.figure(2, figsize=(8, 6))
    for planet in range(1, 5):
        plt.subplot(2, 2, planet)
        plt.plot(np.linspace(0, timestep*num_steps, num=num_steps), normalised_planet_energy[:, planet-1], label=body_names[planet], color=plot_colours[planet])
        plt.legend()
    plt.suptitle("Terrestrial Planets Stability")
    plt.ylabel('Normalised Kinetic Energy')
    plt.xlabel('Time')
    plt.savefig("images/unstructured_3_terrestrial_body_energy_subplot.pdf")

    #gas giants
    plt.figure(3, figsize=(8, 6))
    plt.title("Normalised Gaseous Planet Orbit Stability")
    for planet in range(5, 9):
        plt.subplot(2, 2, planet-4)
        plt.plot(np.linspace(0, timestep*num_steps, num=num_steps), normalised_planet_energy[:, planet-1], label=body_names[planet], color=plot_colours[planet])
        plt.legend()
    plt.ylabel('Normalised Kinetic Energy')
    plt.xlabel('Time')
    plt.legend()
    plt.savefig("images/unstructured_3_gaseous_planet_stability.pdf")
    plt.show()