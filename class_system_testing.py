import math
import numpy as np
import matplotlib.pyplot as plt

AU = 1.495978707e11
G = 6.67e-11
M_SUN = 1.98892e30

#create planet class to make defining new planets easier
class Planet:
    def __init__(self, mass, radius):
        self.mass = mass
        self.radius = radius

    def __str__(self):
        return f'Mass: {self.mass} \nInitial radius: {self.radius}'
        
sun = Planet(1.98892E30, 0) #i know its not a planet shut up
earth = Planet(5.9742E24, 1.745978707E11)
#jupiter = Planet(1.8986E27, 5.2 * 1.495978707e11)

print(earth)

print(earth.mass)
print(earth.radius)

#test here