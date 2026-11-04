import time
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def initial_positions(radii):
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
    a_x_position , a_y_position = vx_array,  vy_array                   #abcd values will all be arrays 
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

#constants
AU = 1.495978707e11
G = 6.67e-11
M_SUN = 1.98892e30
timestep = 3600 #3600 seconds = 1 hour
num_steps = 10 * (365 * 24)  #multiple of earth years

#use arrays to store data for each body
#data taken from https://solarsystem.nasa.gov/planet-compare/ 
body_names = np.array(["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"])
masses = np.array([1.98892e30, 3.30104E23 , 4.86732E24, 5.9742E24, 6.4169E23, 1.89813E27, 5.68319E26, 8.68103E25, 1.0241E26])
initial_radii = np.array([0, 5.7909227E10 , 1.08209475E11, 1.495978707e11, 2.27956E11, 7.78340821E11, 1.426666422E12, 2.870658186E12, 4.498396441E12])

#create blank arrays of equal dimensions to store position and velocity data
x_positions = np.zeros((num_steps,len(body_names))); y_positions = np.zeros((num_steps,len(body_names)))
x_velocities = np.zeros((num_steps,len(body_names))); y_velocities = np.zeros((num_steps,len(body_names)))


#initialise position of all bodies surrounding the sun
x_positions[0], y_positions[0] = initial_positions(initial_radii)

#initialise velocities of all planets ensuring to set sun's velocity to 0
x_velocities[0], y_velocities[0] = initial_velocity(x_positions[0], y_positions[0])

#initialise velocity of sun to give system zero net momentum
x_velocities[0, 0] = -np.sum(masses[1:] * x_velocities[0, 1:]) / masses[0]
y_velocities[0, 0] = -np.sum(masses[1:] * y_velocities[0, 1:]) / masses[0]


#main simulation loop
start_time = time.time()
for loop in tqdm(range(1, num_steps), ncols=90, ascii=" ▖▘▝▗▚▞█", desc="Running Simulation", colour="#00bb82"):  #"░▒█" -- " ▖▘▝▗▚▞█"
    x_positions[loop], y_positions[loop], x_velocities[loop], y_velocities[loop] = body_update(x_positions[loop-1], y_positions[loop-1], x_velocities[loop-1], y_velocities[loop-1], timestep)

#display simulation runtime data to terminal
end_time = time.time()
runtime = end_time - start_time
print(f'Simulation time: {timestep*num_steps} seconds')
print(f'Runtime for {num_steps} iterations: {round(runtime, 5)} seconds')

#plot orbits of bodies to screen
plt.figure(figsize=(6, 6))
for i in range(len(body_names)):
    plt.plot(x_positions[:, i], y_positions[:, i], label=body_names) #legend is not functional
plt.axis('equal') #keeps square aspect ratio
plt.show()