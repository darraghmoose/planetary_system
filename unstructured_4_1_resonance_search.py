import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from orbits import initial_positions, initial_velocity, body_update

#global parameters
AU = 1.495978707e11
G = 6.67e-11
body_names = np.array(["Sun", "Jupiter", "Asteroid"])
masses = np.array([1.98892e30, 1.89813E27, 0.95E21])

def kinetic_energy(body_mass, body_x_velocity, body_y_velocity):
    energy = 0.5 * body_mass * (body_x_velocity**2 + body_y_velocity**2)
    return energy

#main simulation loop
def asteroid_stability_simulation(target_radii):
    """This function accepts an array of radii of asteroids. A simulation is ran for timestep*num_steps seconds which 
    returns relevant information on the positions, velocities, energies and periods of the asteroids input and jupiter's period. """

    #initialise arrays which will be returned to the user
    simulation_asteroid_x_positions = [[] for i in range(len(target_radii))]; simulation_asteroid_y_positions = [[] for i in range(len(target_radii))]
    simulation_asteroid_energies = [[] for i in range(len(target_radii))]; simulation_asteroid_deviations = [[] for i in range(len(target_radii))]
    simulation_asteroid_crossings = [[] for i in range(len(target_radii))]; simulation_jupiter_crossings = [[] for i in range(len(target_radii))]

    for i in range(len(target_radii)):
        initial_radii = np.array([0, 7.78479E11, target_radii[i]])

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
        body_energy = np.zeros((num_steps, len(body_names)))
        body_energy[0] = kinetic_energy(masses[0:], x_velocities[0][0:], y_velocities[0][0:])
        
        crossings = [[] for j in range(len(body_names)-1)] 
        for loop in tqdm(range(1, num_steps), ncols=90, ascii=" ░▒█", desc=f'Asteroid {i+1}/{len(target_radii)}', colour="#00bb82"):  #" ░▒█" -- " ▖▘▝▗▚▞█"
            x_positions[loop], y_positions[loop], x_velocities[loop], y_velocities[loop] = body_update(x_positions[loop-1], y_positions[loop-1], x_velocities[loop-1], y_velocities[:][loop-1], timestep, masses)
            #Measure asteroid KE to gague orbit stability
            body_energy[loop] = kinetic_energy(masses[:len(masses)], x_velocities[loop][:len(masses)], y_velocities[loop][:len(masses)])

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

        #calculate energy deviation
        normalised_energy = body_energy/np.mean(body_energy, axis=0)
        energy_deviation = np.max(np.abs(normalised_energy-1), axis=0)

        #save results for this iteration to asteroid info 
        simulation_asteroid_x_positions[i] = x_positions[:, 2]  #type: ignore
        simulation_asteroid_y_positions[i] = y_positions[:, 2]  #type: ignore
        simulation_asteroid_energies[i] = body_energy[:, 2]     #type: ignore
        simulation_asteroid_deviations[i] = energy_deviation[2]
        simulation_asteroid_crossings[i] = crossings[1]
        simulation_jupiter_crossings[i] = crossings[0]

    return simulation_asteroid_x_positions, simulation_asteroid_y_positions, simulation_asteroid_energies, simulation_asteroid_deviations, simulation_asteroid_crossings, simulation_jupiter_crossings

#simulation parameters
timestep = 24 * 3600    #24 hours is accurate enough for mars and so also asteroid belt.  
num_steps = 50 * (365)  #runtime in years
start_radius = 2.4 *AU
end_radius = 2.6 * AU
radius_increment = 0.01 * AU
num_asteroids = int((end_radius-start_radius)/radius_increment +1)
target_radii = np.linspace(start_radius, end_radius, num_asteroids) #if want to simulate only certain radii, replace target_radii with list of desired radii


#run simulation
asteroid_x_positions, asteroid_y_positions, asteroid_energies, asteroid_deviations, asteroid_crossings, jupiter_crossings = asteroid_stability_simulation(target_radii)


#display results 
print(f'Orbital period results:\n')
print(f'{'Radius (AU)':<15}{'Seconds':>15}{'Days':>15}{'Resonance':>15}{'Deviation (%)':>15}')
print(f'---------------------------------------------------------------------------')
for i in range(len(target_radii)):
    asteroid_times = asteroid_crossings[i] #number of times the asteroid orbits the sun
    jupiter_times = jupiter_crossings[i]
    if len(asteroid_times) < 2 or len(jupiter_times) < 2: #body needs to pass x axis twice to measure period
        print(f'{f'{i}':<15}{'Insufficient crossings':>75}')
        continue
    orbital_period = np.diff(asteroid_times).mean()
    jupiter_period = np.diff(jupiter_times).mean()
    period_ratio = jupiter_period / orbital_period
    print(f'{target_radii[i]/AU:<15,.2f}{orbital_period:>15,.1f}{orbital_period/3600:>15,.2f}{period_ratio:>15,.3f}{asteroid_deviations[i]*100:15,.1f}')
