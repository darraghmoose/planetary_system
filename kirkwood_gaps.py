import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def body_initial_positions(radii):
    #intake radius of planets orbit and return random x,y coords to start planet at
    theta = np.random.uniform(0, 2*np.pi, size = len(radii))    #theta in radians
    random_x = radii * np.cos(theta)   
    random_y = radii * np.sin(theta)
    #ensure sun remains at (0,0)
    random_x[0], random_y[0] = (0.0, 0.0)

    return random_x, random_y

def initial_velocity(x_start: np.ndarray, y_start: np.ndarray):
    #returns initial velocity of planets at the start of simulaiton based only on orbital radius
    initial_x_velocity = np.zeros(len(x_start))
    initial_y_velocity = np.zeros(len(y_start))
    orbit_radius = np.sqrt(x_start**2 + y_start**2)

    for i in range(len(orbit_radius)):
        if orbit_radius[i] == 0:
            continue
        speed = np.sqrt(G * masses[0] / orbit_radius[i])
        initial_x_velocity[i] = speed * y_start[i] / orbit_radius[i]
        initial_y_velocity[i] = -speed * x_start[i] / orbit_radius[i]
    
    return initial_x_velocity, initial_y_velocity

def acceleration(x_positions, y_positions) -> tuple:
    """ Function accepts x and y positions of all bodies at once and calculates the acceleration applied to all bodies simultaneously. 
    Returns updated acceleration for all bodies"""
    x_seperations = x_positions[None, :] - x_positions[:, None] #create 2D array where item[i][j] is the distance between j and i (xj - xi)
    y_seperations = y_positions[None, :] - y_positions[:, None] #same array but for y position
    #"https://stackabuse.com/guide-to-numpy-matrix-subtraction/" info on how vector subreaction works to make difference matrix
    seperations = np.sqrt(x_seperations**2 + y_seperations**2) #calculate geometric seperation between bodies. works itemwise across both arrays
    """when i=j their seperation =0. this iwll lead to 0/0 error when calculating k. We want k=0 for this case. 
    To do this, we use fill_diagonal() to set seperation[i][i]=infinity to give 0/infinity which the computer will read as 0 giving k=0 """
    np.fill_diagonal(seperations, np.inf) #np.inf gives infinity to guarantee k=0

    #calculate new accelerations of all bodies and save to 1D array
    x_acceleration = (G * np.sum(masses[None, :] * x_seperations / (seperations**3), axis=1)) #np.sum ensures the acceleration change is imparted on the body by all other bodies on the system at the same time
    y_acceleration = (G * np.sum(masses[None, :] * y_seperations / (seperations**3), axis=1)) #axis=1 ensures dimensions are right for returning as 1D array
    return x_acceleration, y_acceleration

def body_update(x_array, y_array, vx_array, vy_array, dt: float) -> tuple:
    """
    Accepts arrays of positions and velocities
    Returns updated position and velocity of subject in x and y
    """

    #calculate all a values
    a_x_position , a_y_position = vx_array,  vy_array             #abcd values will all be arrays 
    a_x_velocity, a_y_velocity = acceleration(x_array, y_array)   #updated velocity is an acceleration by definition of runge kutta method, that is reason function is named acceleration

    #calculate all b values
    b_x_position, b_y_position = vx_array + 0.5*dt*a_x_velocity, vy_array + 0.5*dt*a_y_velocity
    b_x_velocity, b_y_velocity = acceleration(x_array + 0.5*dt*a_x_position, y_array + 0.5*dt*a_y_position)

    #calculate all c values
    c_x_position, c_y_position = vx_array + 0.5*dt*b_x_velocity, vy_array + 0.5*dt*b_y_velocity
    c_x_velocity, c_y_velocity = acceleration(x_array + 0.5*dt*b_x_position, y_array + 0.5*dt*b_y_position)
   
    #calculate all d values
    d_x_position, d_y_position = vx_array + dt*c_x_velocity, vy_array + dt*c_y_velocity
    d_x_velocity, d_y_velocity = acceleration(x_array + dt*c_x_position, y_array + dt*c_y_position)

    #calculate weighted average of values to update all variables
    updated_x_positions = x_array + ((dt/6) * (a_x_position + 2*b_x_position + 2*c_x_position + d_x_position))
    updated_y_positions = y_array + ((dt/6) * (a_y_position + 2*b_y_position + 2*c_y_position + d_y_position))
    updated_x_velocities = vx_array + ((dt/6) * (a_x_velocity + 2*b_x_velocity + 2*c_x_velocity + d_x_velocity))
    updated_y_velocities = vy_array + ((dt/6) * (a_y_velocity + 2*b_y_velocity + 2*c_y_velocity + d_y_velocity))
    return updated_x_positions, updated_y_positions, updated_x_velocities, updated_y_velocities

def kinetic_energy(body_mass, body_x_velocity, body_y_velocity):
    energy = 0.5 * body_mass * (body_x_velocity**2 + body_y_velocity**2)
    return energy

#constants
AU = 1.495978707e11
G = 6.67e-11
M_SUN = 1.98892e30

#simulation runtime parameters
timestep = 24 * 3600        #seconds per timestep
num_steps = 4 * (365 * 24)  #iterations per simulaiton

#use arrays to store data for each body      -      data taken from https://solarsystem.nasa.gov/planet-compare/ 
"""
body_names = np.array(["Sun", "Jupiter", "5:1 Resonance", "4:1 Resonance", "3:1 Resonance", "Asteroid Belt"])
masses = np.array([1.98892e30, 1.89813E27, 0.95E21, 0.95E21, 95E21, 95E21])
initial_radii = np.array([0, 5.204*AU, 1.780*AU, 2.065*AU, 2.502*AU, 2.8*AU])

body_names = np.array(["Sun", "Jupiter", "Unstable 1", "Unstable 2", "Asteroid Belt"])
masses = np.array([1.98892e30, 1.89813E27, 0.95E21, 0.95E21, 95E21])
initial_radii = np.array([0, 5.2*AU, (3.28 * AU), (4 * AU), 2.8*AU])       #initial_radii = np.array([0, 5.2*AU, (3.276 * AU), (4 * AU), 2.8*AU])


body_names = np.array(["Sun", "Jupiter", "Launched"])
masses = np.array([1.98892e30, 1.89813E27, 0.95E21])
initial_radii = np.array([0, 5.2*AU, 4.5*AU])
"""

#Gravity assist
body_names = np.array(["Sun", "Jupiter", "Asteroid 1", "Asteroid 2", "Asteroid 3", "Asteroid 4", "Asteroid 5"])
masses = np.array([1.98892e30, 1.89813E27, 0.95E21, 0.95E21, 95E21, 95E21, 95E21])
initial_radii = np.array([0, 5.2*AU, 3.96* AU, 3.98* AU, 4.00 * AU, 4.02*AU, 4.04*AU])  


plot_colours = ['#ef0000', '#c98226', '#8c8987', '#dfb18c', '#3c6690', '#ff6b4d', '#f8d280', '#5d8b8e', '#4e5a74']

#create blank arrays of equal dimensions to store position and velocity data
x_positions = np.zeros((num_steps,len(body_names))); y_positions = np.zeros((num_steps,len(body_names)))
x_velocities = np.zeros((num_steps,len(body_names))); y_velocities = np.zeros((num_steps,len(body_names)))

#initialise position and velocity of all bodies surrounding the sun
x_positions[0], y_positions[0] = body_initial_positions(initial_radii)
x_velocities[0], y_velocities[0] = initial_velocity(x_positions[0], y_positions[0])

#initialise velocity of sun to give system zero net momentum
x_velocities[0, 0] = -np.sum(masses[1:] * x_velocities[0, 1:]) / masses[0]
y_velocities[0, 0] = -np.sum(masses[1:] * y_velocities[0, 1:]) / masses[0]

#track KE of asteroids
asteroid_energy = np.zeros((num_steps, len(body_names)-2))
asteroid_energy[0] = kinetic_energy(masses[2:], x_velocities[0][2:], y_velocities[0][2:])

#main simulation loop
crossings = [[] for i in range(len(body_names)-1)]
for loop in tqdm(range(1, num_steps), ncols=90, ascii=" ░▒█", desc="Asteroid Simulation", colour="#00bb82"):  #" ░▒█" -- " ▖▘▝▗▚▞█"
    x_positions[loop], y_positions[loop], x_velocities[loop], y_velocities[loop] = body_update(x_positions[loop-1], y_positions[loop-1], x_velocities[loop-1], y_velocities[:][loop-1], timestep)
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
deviation = np.max(normalised_energy, axis=0)# - np.min(normalised_energy, axis=0) 
print(deviation)

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
plt.figure(1, figsize=(6, 6))
for i in range(len(body_names)):
    plt.plot(x_positions[:, i], y_positions[:, i], color=plot_colours[i], label=body_names[i])
plt.xlabel('X [m]')
plt.ylabel('Y [m]')
plt.legend()

#create normalised plot of asteroid energy
plt.figure(2, figsize=(6, 6))
plt.title("Normalised Asteroid Stability")
for asteroid in range(len(normalised_energy[0])):
    plt.plot(np.linspace(0, timestep*num_steps, num=num_steps), normalised_energy[:, asteroid], color=plot_colours[asteroid+2], label=body_names[asteroid+2]) 
plt.ylabel('Normalised Kinetic Energy')
plt.xlabel('Time')
plt.legend()

#display plots
plt.show()
"""
TODO
Rewrite kirkwood gaps programme to iterate through individual astroidsto avoid nearby asteroid affecting each other
This can be used to determine radii and resonance periods for these orbits
"""