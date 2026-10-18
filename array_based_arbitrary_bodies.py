import math
import numpy as np
import time
import matplotlib.pyplot as plt

def acceleration(subject_x: float, subject_y: float, affecting_x: float, affecting_y: float, affecting_index: int) -> tuple:
    
    #calculate distance between both bodies
    seperation = math.sqrt((subject_x - affecting_x)**2 + (subject_y - affecting_y)**2)
    #calculate constant to multiply (xj-xi)     
    k = -G * masses[affecting_index] / seperation**3

    return k * (affecting_x - subject_x), k * (affecting_y - subject_y) * k

def body_update(subject_index: int, affecting_index: int, timestep: float) -> tuple:
    """
    Accepts the index of subject body whose velocity and position is being updated, index of affecting body who is pulling on subject and the timestep
    Returns updated position and velocity of subject in x and y

    modififcation:
    accept indicies of bodies in question - subject body and affecting body
    calculate how affecting body changes position and velocity of subject body
    return updated position and velocity of subject body
    """
    #find position and velocity of subject body and position of affecting body
    #Subject body
    x_position, y_position = x_positions[subject_index], y_positions[subject_index]
    x_velocity, y_velocity = x_velocities[subject_index], y_velocities[subject_index]

    #Position of body affecting subject body
    affecting_body_x_position, affecting_body_y_position = x_positions[affecting_index], y_positions[affecting_index]


    #calculate all a values
    a_x_position , a_y_position = x_velocity,  y_velocity           #updates both positions at once
    a_x_velocity, a_y_velocity = acceleration(x_position, y_position, affecting_body_x_position, affecting_body_y_position, affecting_index)

    #calculate all b values
    b_x_position, b_y_position = x_velocity + 0.5*timestep*a_x_velocity, y_velocity + 0.5*timestep*a_y_velocity
    b_x_velocity, b_y_velocity = acceleration(x_position + 0.5*timestep*a_x_position, y_position + 0.5*timestep*a_y_position, affecting_body_x_position, affecting_body_y_position, affecting_index)

    #calculate all c values
    c_x_position, c_y_position = x_velocity + 0.5*timestep*b_x_velocity, y_velocity + 0.5*timestep*b_y_velocity
    c_x_velocity, c_y_velocity = acceleration(x_position + 0.5*timestep*b_x_position, y_position + 0.5*timestep*b_y_position, affecting_body_x_position, affecting_body_y_position, affecting_index)
   
    #calculate all d values
    d_x_position, d_y_position = x_velocity + timestep*c_x_velocity, y_velocity + timestep*c_y_velocity
    d_x_velocity, d_y_velocity = acceleration(x_position + timestep*c_x_position, y_position + timestep*c_y_position, affecting_body_x_position, affecting_body_y_position, affecting_index)

    #calculate weighted average of values to update all variables
    updated_x_position = x_position + ((timestep/6) * (a_x_position + 2*b_x_position + 2*c_x_position + d_x_position))
    updated_y_position = y_position + ((timestep/6) * (a_y_position + 2*b_y_position + 2*c_y_position + d_y_position))
    updated_x_velocity = x_velocity + ((timestep/6) * (a_x_velocity + 2*b_x_velocity + 2*c_x_velocity + d_x_velocity))
    updated_y_velocity = y_velocity + ((timestep/6) * (a_y_velocity + 2*b_y_velocity + 2*c_y_velocity + d_y_velocity))
    return updated_x_position, updated_y_position, updated_x_velocity, updated_y_velocity

#useful constants
AU = 1.495978707e11
G = 6.67e-11
M_SUN = 1.98892e30
timestep = 100
simulation_time = timestep * 10000

#use arrays to store data for each body
body_names = ["Sun", "Earth"]
x_positions = np.array([0, 1*AU])
y_positions = np.array([0, 0])
x_velocities = np.array([0, 0])
y_velocities = np.array([0, -1 * math.sqrt(G * M_SUN / AU )])
masses = np.array([1.98892e30, 5.9742E24])

#create array to store data which can be plotted
positions_plotting = np.zeros((len(body_names), simulation_time))   #need to update to store x and y data

x_testing_array = [1*AU]
y_testing_array = [0]

#add initial positions to positions_plotting array



#main loop to run simulation
start_time = time.time()
for loop in range(simulation_time):
    #below loop does 1 iteration per entire loop
    for i in range(len(body_names)):        #how does body i move
        for j in range(len(body_names)):    #from the efefcts of body j
            if j != i:                      #it cannot affect itself
                #run update body but need to modify function first
                x_positions[i], y_positions[i], x_velocities[i], y_velocities[i] = body_update(i, j, timestep=timestep)
                if i == 1:
                    x_testing_array.append(x_positions[i])
                    y_testing_array.append(y_positions[i])



#timing testing
end_time = time.time()
runtime = end_time - start_time
print(f'Runtime for {simulation_time} iterations: {round(runtime, 4)} seconds')
print(f'Average simulation step time: {runtime/simulation_time}')


plt.plot(x_testing_array, y_testing_array)
plt.show()
#dont like the look of the plot


"""
Plotting:
need to save position data for all bodies across time to plot
consider using 2D array 
#columns = #planets
#rows = iteration at each time

gives array of size(#planets, simulation_time)
will need to initialise with positions at time=0

"""