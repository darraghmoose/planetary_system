import numpy as np
from tqdm import tqdm
from orbits import initial_positions, initial_velocity, body_update

#simulation settings
num_steps = int(3 * (365 * 24))

#Orbiting bodies data
body_names = np.array(["Sun", "Mars"])
masses = np.array([1.98892e30, 5.9742E24])
initial_radii = np.array([0, 2.27956E11])

#main simulation function
def run_simulation(timestep):
    #create blank arrays of equal dimensions to store position and velocity data
    x_positions = np.zeros((num_steps,len(body_names))); y_positions = np.zeros((num_steps,len(body_names)))
    x_velocities = np.zeros((num_steps,len(body_names))); y_velocities = np.zeros((num_steps,len(body_names)))

    #initialise position of all bodies surrounding the sun
    x_positions[0], y_positions[0] = initial_positions(initial_radii)

    #initialise velocities of all planets
    x_velocities[0], y_velocities[0] = initial_velocity(x_positions[0], y_positions[0], masses)

    #main simulation loop
    crossings = [[] for i in range(len(body_names)-1)]
    for loop in range(1, num_steps):
        x_positions[loop], y_positions[loop], x_velocities[loop], y_velocities[loop] = body_update(x_positions[loop-1], y_positions[loop-1], x_velocities[loop-1], y_velocities[:][loop-1], timestep, masses)
        #check for axis crossing
        for planet in range(1, len(body_names)):
            y_coord_now = y_positions[loop][planet] - y_positions[loop][0]
            y_coord_prev = y_positions[loop-1][planet] - y_positions[loop-1][0]
            x_coord_now = x_positions[loop][planet] - x_positions[loop][0]

            #look for sign change across axis 
            if y_coord_prev * y_coord_now < 0 and x_coord_now > 0:
                #timestep at [loop] = loop * timesteps
                time_after = loop * timestep
                crossing_time = time_after - (y_coord_now)*(timestep / (y_coord_now - y_coord_prev))
                crossings[planet-1].append(crossing_time)

    #calculate period from crossings
    for i in range(1, len(body_names)):
        times = crossings[i-1]
        orbital_period = np.diff(times).mean()

    return orbital_period

#run simulation on range of timesteps
timestep_array = [3600*1, 3600*2, 3600*4,3600*6, 3600*12, 3600*24, 3600*72, 3600*24*7, 3600*24*14]
periods_array = []
for i in tqdm(range(len(timestep_array)), ncols=90, ascii=" ░▒█", desc="Running Simulation", colour="#00bb82"):
    periods_array.append(run_simulation(timestep_array[i]))

#display results
print(f'Orbital period results:\n')
print(f'{'Timestep (hours)':<18}{'Seconds':>15}{'Hours':>15}{'Days':>15}')
print(f'---------------------------------------------------------')
for i in range(len(timestep_array)):
    print(f'{timestep_array[i]/3600:<18,}{periods_array[i]:>15,.1f}{periods_array[i]/3600:>15,.2f}{periods_array[i]/(3600*24):>15,.2f}')