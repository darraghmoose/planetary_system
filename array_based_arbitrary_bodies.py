import math
import numpy as np
import time
import matplotlib.pyplot as plt
#from tqdm import tqdm

def acceleration(subject_x: float, subject_y: float, affecting_x: float, affecting_y: float, affecting_index: int) -> tuple:
    
    #calculate distance between both bodies
    seperation = math.sqrt((subject_x - affecting_x)**2 + (subject_y - affecting_y)**2)
    #calculate constant to multiply (xj-xi)     
    k = -G * masses[affecting_index] / seperation**3

    return k * (affecting_x - subject_x), k * (affecting_y - subject_y) 

def body_update(subject_index: int, affecting_index: int, loop: int , timestep: float) -> tuple:
    """
    Accepts the index of subject body whose velocity and position is being updated, index of affecting body who is pulling on subject and the timestep
    Returns updated position and velocity of subject in x and y
    """
    #find position and velocity of subject body and position of affecting body
    #Subject body
    prev_loop = loop - 1
    x_position, y_position = x_positions[prev_loop][subject_index], y_positions[prev_loop][subject_index]
    x_velocity, y_velocity = x_velocities[prev_loop][subject_index], y_velocities[prev_loop][subject_index]

    #Position of body affecting subject body
    affecting_body_x_position, affecting_body_y_position = x_positions[prev_loop][affecting_index], y_positions[prev_loop][affecting_index]


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
num_steps = 6
simulation_time = timestep * num_steps #time at simulation step i = i*timestep

#use arrays to store data for each body
body_names = ["Sun", "Earth"]
masses = np.array([1.98892e30, 5.9742E24])
#create arrays of equal dimensions to store position and velocity data
x_positions = np.zeros((num_steps,len(body_names)))
y_positions = np.zeros((num_steps,len(body_names)))
x_velocities = np.zeros((num_steps,len(body_names)))
y_velocities = np.zeros((num_steps,len(body_names)))
#initialise these arrays with system's starting positions
x_positions[0] = np.array([0, 1*AU]); print(f'X-position[0][1] = {x_positions[0][1]}')
y_positions[0] = np.array([0, 0])
x_velocities[0] = np.array([0, 0])
y_velocities[0] = np.array([0, -1 * math.sqrt(G * M_SUN / AU )])


#testing
x_testing_array = [1*AU]
y_testing_array = [0]

#main loop to run simulation
start_time = time.time()
for loop in range(1, num_steps): #begin at index [1] as initial positions have already been declared
    #below loop does 1 iteration per entire loop
    for i in range(len(body_names)):        #how does body i move
        for j in range(1):#len(body_names)):    #from the efefcts of body j
            if j != i:                      #it cannot affect itself
                #run update body but need to modify function first
                next_x_position, next_y_position, next_x_velocity, next_y_velocity = body_update(i, j, loop, timestep=timestep)
                #save next positions to array structure
                x_positions[loop][i] = next_x_position
                y_positions[loop][i] = next_y_position                
                x_velocities[loop][i] = next_x_velocity             
                y_velocities[loop][i] = next_y_velocity


#timing testing
end_time = time.time()
runtime = end_time - start_time
print(f'Runtime for {simulation_time/timestep} iterations: {round(runtime, 4)} seconds')

print(x_positions.shape)

for i in range(0, num_steps):
    print(x_positions[i,1])



#print(len(x_positions[0:num_steps, 1]))
plt.plot(x_positions[0:num_steps, 1], y_positions[0:num_steps, 1])
plt.show()
#plot wrong </3

"""
Plotting:
need to save position data for all bodies across time to plot
consider using 2D array 
#columns = #planets
#rows = iteration at each time

gives array of size(#planets, simulation_time)
will need to initialise with positions at time=0

"""