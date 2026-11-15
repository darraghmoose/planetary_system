import numpy as np

G = 6.67e-11
M_SUN=1.98892e30

def initial_positions(radii: np.ndarray) -> tuple:
    #intake radius of planets orbit and return random x,y coords to start planet at
    theta = np.random.uniform(0, 2*np.pi, size = len(radii))    #theta in radians
    random_x = radii * np.cos(theta)   
    random_y = radii * np.sin(theta)
    #ensure sun remains at (0,0)
    random_x[0], random_y[0] = (0.0, 0.0)

    return random_x, random_y

def initial_velocity(x_start: np.ndarray, y_start: np.ndarray, mass):
    #returns initial velocity of planets at the start of simulaiton based only on orbital radius
    initial_x_velocity = np.zeros(len(x_start))
    initial_y_velocity = np.zeros(len(y_start))
    orbit_radius = np.sqrt(x_start**2 + y_start**2)

    for i in range(len(orbit_radius)):
        if orbit_radius[i] == 0: #ignore the sun initially, it is given velocity to give system momentum = 0
            continue
        speed = np.sqrt(G * mass[0] / orbit_radius[i])
        initial_x_velocity[i] = speed * y_start[i] / orbit_radius[i]
        initial_y_velocity[i] = -speed * x_start[i] / orbit_radius[i]
    
    return initial_x_velocity, initial_y_velocity

def elliptical_initial_velocity(x_start, y_start, eccentricity):
    """Accepts the starting position of each body at apogee and its eccentricity. 
    Returns the initial velocity needed for that orbit."""
    sun_mass = 1.98892e30 #kg
    initial_x_velocity = np.zeros(len(x_start))
    initial_y_velocity = np.zeros(len(y_start))
    apogee_radius = np.hypot(x_start, y_start)

    for i in range(len(apogee_radius)):
        if apogee_radius[i] == 0.0: #ignore the sun initially, it is later given velocity to give system momentum = 0
            continue
        speed = np.sqrt(G * sun_mass * (1 - eccentricity[i]) / apogee_radius[i])
        initial_x_velocity[i] = speed * y_start[i] / apogee_radius[i]
        initial_y_velocity[i] = -speed * x_start[i] / apogee_radius[i]

    return initial_x_velocity, initial_y_velocity

def orbital_elements(x, y, vx, vy, xs, ys, vxs, vys, mu=G*M_SUN):
    """Returns semi-major axis and eccentricity of body relative to the sun."""
    rx, ry = x - xs, y - ys
    rvx, rvy = vx - vxs, vy - vys
    r = np.hypot(rx, ry)
    v2 = rvx**2 + rvy**2
    a = 1 / (2/r - v2/mu)
    h = rx*rvy - ry*rvx     #specific angular momentum
    e = np.sqrt(np.maximum(0, 1 + (v2 - 2*mu/r) * h**2 / mu**2))
    return a, e

def acceleration(x_positions, y_positions, mass) -> tuple:
    """ Function accepts x and y positions of all bodies at once and calculates the acceleration applied to all bodies simultaneously. 
    Returns updated acceleration for all bodies"""
    x_seperations = x_positions[None, :] - x_positions[:, None] #create 2D array where item[i][j] is the distance between j and i (xj - xi)
    y_seperations = y_positions[None, :] - y_positions[:, None] #same array but for y position
    #"https://stackabuse.com/guide-to-numpy-matrix-subtraction/" info on how vector subreaction works to make difference matrix
    seperations = np.sqrt(x_seperations**2 + y_seperations**2) #calculate seperation between bodies
    """when i=j their seperation =0. this will lead to 0/0 error when calculating k. We want k=0 for this case. 
    To do this, we use fill_diagonal() to set seperation[i][i]=infinity to give 0/infinity which the computer will read as 0 giving k=0 """
    np.fill_diagonal(seperations, np.inf) #np.inf gives infinity to guarantee k=0

    #calculate new accelerations of all bodies and save to 1D array
    x_acceleration = (G * np.sum(mass[None, :] * x_seperations / (seperations**3), axis=1)) #np.sum ensures the acceleration change is imparted on the body by all other bodies on the system at the same time
    y_acceleration = (G * np.sum(mass[None, :] * y_seperations / (seperations**3), axis=1)) #axis=1 ensures dimensions are right for returning as 1D array
    return x_acceleration, y_acceleration

def body_update(x, y, vx, vy, dt: float, mass) -> tuple:
    """
    Accepts arrays of positions and velocities
    Returns updated position and velocity of subject in x and y
    """

    #calculate all a values
    a_x_position , a_y_position = vx,  vy             #abcd values will all be arrays 
    a_x_velocity, a_y_velocity = acceleration(x, y, mass)   #updated velocity is an acceleration by definition of runge kutta method, that is reason function is named acceleration

    #calculate all b values
    b_x_position, b_y_position = vx + 0.5*dt*a_x_velocity, vy + 0.5*dt*a_y_velocity
    b_x_velocity, b_y_velocity = acceleration(x + 0.5*dt*a_x_position, y + 0.5*dt*a_y_position, mass)

    #calculate all c values
    c_x_position, c_y_position = vx + 0.5*dt*b_x_velocity, vy + 0.5*dt*b_y_velocity
    c_x_velocity, c_y_velocity = acceleration(x + 0.5*dt*b_x_position, y + 0.5*dt*b_y_position, mass)
   
    #calculate all d values
    d_x_position, d_y_position = vx + dt*c_x_velocity, vy + dt*c_y_velocity
    d_x_velocity, d_y_velocity = acceleration(x + dt*c_x_position, y + dt*c_y_position, mass)

    #calculate weighted average of values to update all variables
    updated_x_positions = x + ((dt/6) * (a_x_position + 2*b_x_position + 2*c_x_position + d_x_position))
    updated_y_positions = y + ((dt/6) * (a_y_position + 2*b_y_position + 2*c_y_position + d_y_position))
    updated_x_velocities = vx + ((dt/6) * (a_x_velocity + 2*b_x_velocity + 2*c_x_velocity + d_x_velocity))
    updated_y_velocities = vy + ((dt/6) * (a_y_velocity + 2*b_y_velocity + 2*c_y_velocity + d_y_velocity))

    return updated_x_positions, updated_y_positions, updated_x_velocities, updated_y_velocities