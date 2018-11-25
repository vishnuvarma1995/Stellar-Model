# Static Stellar Model
import numpy as np
from matplotlib import pyplot as plt
# Starting equations:
    # Mass continuity equations: dm/dr = 4*pi*r**2*roh
    # Hydrostatic Equilibrium: dP/dr = roh*G*m/R**2
# Major assumptions:
    # The entire star is an ideal gas
    # Density varies linearly
# Constants
roh_c = 162200 #kg/m^3
Star_Radius = 695700000 #m
G = 6.67408*(10**-11) #m^3 kg-1 s-2
mu = 1/2
mass_H = 1.6737236*(10**-27) #kg
boltzmann = 1.38064852*(10**-23)# m2 kg s-2 K-1
radius = []
mass = []
density = []
temperature = []
pressure = []
time =[]
Total_Mass = (np.pi/3)*roh_c*Star_Radius**3
Central_Pressure = (5*G*Total_Mass**2)/(4*np.pi*Star_Radius**4)
# Central_Density = (3*Total_Mass/np.pi)/(Star_Radius**3)

# Iterate to get list of radius values
r = Star_Radius/1000
for i in range(0, Star_Radius, 10000):
    radius.append(i)

def Mass(radius_list):
    const = 4*np.pi*roh_c
    for r in radius_list:
        m = const*((r**3)/3 - (r**4)/(4*Star_Radius))
        mass.append(m)

def Density(radius_list):
    for r in radius_list:
        roh = roh_c*(1-r/Star_Radius)  
        density.append(roh)

def Pressure(radius_list):
    const = 4*np.pi*G*(roh_c**2)
    for r in radius_list:
        r_dependence = (r**2)/6 - (7*r**3)/(36*Star_Radius) + (r**4)/(16*Star_Radius**2)
        P = Central_Pressure - const*r_dependence
        pressure.append(P)

def Temperature(density_list, pressure_list):
    const = (mu*mass_H)/boltzmann
    for i in range(len(density_list)):
        T = const*(1/density_list[i])*pressure_list[i]
        temperature.append(T)


Mass(radius)
Density(radius)
Pressure(radius)
Temperature(density, pressure)


# Plotting
plt.figure(1) # Mass v Radius
plt.title("Mass against radius")
plt.xlabel("Radius[m]")
plt.ylabel("Mass[kg]")
plt.plot(radius, mass, label='mass', color='red')

plt.figure(2) # Density v Radius
plt.title("Density against radius")
plt.xlabel("Radius[m]")
plt.ylabel("Density[kg/m^3]")
plt.plot(radius, density, label='density', color='blue')

plt.figure(3) # Pressure v Radius
plt.title("Pressure against radius")
plt.xlabel("Radius[m]")
plt.ylabel("Pressure")
plt.plot(radius, pressure, label='pressure', color='green')

plt.figure(4) # Temperature v Radius
plt.title("Temperature against radius")
plt.xlabel("Radius[m]")
plt.ylabel("Temperature[K]")
plt.plot(radius, temperature, label='temperature', color='k')

plt.show()