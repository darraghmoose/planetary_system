import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from orbits import initial_positions, initial_velocity, body_update

#Orbiting bodies data
AU = 1.495978707e11
body_names = np.array(["Sun", "Erid"])
masses = np.array([1.98892e30, 5.9742E24]) #arbitrary planet is given earth mass
initial_radii = np.array([0, (0.4 + 1054/25000) * AU])

def run_simulation(timestep):
    #simulation settings
    num_steps = int(3 * (365 * 24)) 

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

def unstructured3_1():
    #constants & settings
    AU = 1.495978707e11
    timestep = 3600 #3600 seconds = 1 hour
    num_steps = 3 * (365 * 24)  #multiple of earth years

    #create blank arrays of equal dimensions to store position and velocity data
    x_positions = np.zeros((num_steps,len(body_names))); y_positions = np.zeros((num_steps,len(body_names)))
    x_velocities = np.zeros((num_steps,len(body_names))); y_velocities = np.zeros((num_steps,len(body_names)))

    #initialise position of all bodies surrounding the sun
    x_positions[0], y_positions[0] = initial_positions(initial_radii)

    #initialise velocities of all planets
    x_velocities[0], y_velocities[0] = initial_velocity(x_positions[0], y_positions[0], masses)


    #main simulation loop
    for loop in tqdm(range(1, num_steps), ncols=90, ascii=" ░▒█", desc="Plotting orbit", colour="#00bb82"):  #"░▒█" -- " ▖▘▝▗▚▞█"
        x_positions[loop], y_positions[loop], x_velocities[loop], y_velocities[loop] = body_update(x_positions[loop-1], y_positions[loop-1], x_velocities[loop-1], y_velocities[:][loop-1], timestep, masses)


    #plot orbits of bodies to screen
    plt.figure(figsize=(6, 6))
    plot_colours = ['#ef0000', "#3b8eb4"]
    for i in range(len(body_names)):
        if i == 0: #need to use marker for sun as it never moves
            plt.plot(x_positions[0, i]/AU, y_positions[0, i]/AU, 'o', color=plot_colours[i], markersize=12, label=body_names[i]) #[0, i] used to avoid drawing same point over and over
        else:
            plt.plot(x_positions[:, i]/AU, y_positions[:, i]/AU, color=plot_colours[i], label=body_names[i])
    plt.axis('equal')
    plt.xlabel('X TESTING(AU)')
    plt.ylabel('Y (AU)')
    plt.legend()
    plt.savefig("images/structured3.pdf")
    #plt.show()

def unstructured3_2():
    #run simulation on range of timesteps
    timestep_array = [3600*1, 3600*2, 3600*4,3600*6, 3600*12, 3600*24, 3600*72, 3600*24*7, 3600*24*14]
    periods_array = []
    for i in tqdm(range(len(timestep_array)), ncols=90, ascii=" ░▒█", desc="Calculating convergence", colour="#00bb82"):
        periods_array.append(run_simulation(timestep_array[i]))

    #display results
    print(f'Orbital period results:\n')
    print(f'{'Timestep (hours)':<18}{'Seconds':>15}{'Hours':>15}{'Days':>15}')
    print(f'---------------------------------------------------------')
    for i in range(len(timestep_array)):
        print(f'{timestep_array[i]/3600:<18,}{periods_array[i]:>15,.1f}{periods_array[i]/3600:>15,.2f}{periods_array[i]/(3600*24):>15,.2f}')

#run simulation
if __name__=="__main__":
    unstructured3_1()
    unstructured3_2()

"""
TODO
Tidy run_simulation() to be passed into unstructured2_1() also
"""