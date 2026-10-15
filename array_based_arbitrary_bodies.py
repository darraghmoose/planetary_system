import math
import numpy as np

#useful constants
AU = 1.495978707e11
G = 6.67e-11
M_SUN = 1.98892e30
timestep = 100

#use arrays to store data for each body
body_names = ["Sun", "Earth"]
x_positions = np.array([0, 0])
y_positions = np.array([0, 0])
x_velocities = np.array([0, 0])
y_velocities = np.array([0, 0])
masses = np.array([1.98892e30, 5.9742E24])

def update_body(subject_index, x_position_array, y_position_array, x_velocity_array, y_velocity_array, timestep):
    """
    this function is not yet functional.
    can set an r value at the beginning to reuse or does this need to be reupdated with a,b,c,d values as function progresses?
    """
    seperation = math.sqrt(x_position_array**2 + y_position_array**2)

    #calculate all a values
    a_x_position = x_velocity_array 
    a_y_position = y_velocity_array 
    a_x_velocity = (-1 * G * M_SUN * x_position_array) / (seperation)**3                                           
    a_y_velocity = (-1 * G * M_SUN * y_position_array) / (seperation)**3  

    #calculate all b values
    b_x_position = x_velocity_array 
    b_y_position = y_velocity_array 
    b_x_velocity = (-1 * G * M_SUN * (x_position_array + 0.5 * timestep * a_x_velocity)) / (seperation)**3    
    b_y_velocity = (-1 * G * M_SUN * (y_position_array + 0.5 * timestep * a_y_velocity)) / (seperation)**3      
  
    #calculate all c values
    c_x_position = x_velocity_array 
    c_y_position = y_velocity_array 
    c_x_velocity = (-1 * G * M_SUN * (x_position_array + 0.5 * timestep * b_x_velocity)) / (seperation)**3   
    c_y_velocity = (-1 * G * M_SUN * (y_position_array + 0.5 * timestep * b_y_velocity)) / (seperation)**3      
   
    #calculate all d values
    d_x_position = x_velocity_array 
    d_y_position = y_velocity_array 
    d_x_velocity = (-1 * G * M_SUN * (x_position_array + timestep * c_x_velocity)) / (seperation)**3 
    d_y_velocity = (-1 * G * M_SUN * (y_position_array + timestep * c_y_velocity)) / (seperation)**3

    #calculate weighted average of values to update all variables
    updated_x_position = x_position_array + ((timestep/6) * (a_x_position + 2*b_x_position + 2*c_x_position + d_x_position))
    updated_y_position = y_position_array + ((timestep/6) * (a_y_position + 2*b_y_position + 2*c_y_position + d_y_position))
    updated_x_velocity = x_velocity_array + ((timestep/6) * (a_x_velocity + 2*b_x_velocity + 2*c_x_velocity + d_x_velocity))
    updated_y_velocity = y_velocity_array + ((timestep/6) * (a_y_velocity + 2*b_y_velocity + 2*c_y_velocity + d_y_velocity))
    return updated_x_position, updated_y_position, updated_x_velocity, updated_y_velocity