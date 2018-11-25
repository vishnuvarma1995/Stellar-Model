import numpy as np
import matplotlib.pyplot as plt
# Constants
G=6.67E-11
k_B=1.38E-23
m_H=1.67E-27
mu=0.62
Pi=3.141592
# Solar values
Msun=1.989E30
Rsun=7.E8
Rstar=1
Mstar=1
Rstar=Rstar*Rsun # Convert radius to solar radius
Mstar=Mstar*Msun # Convert mass to solar mass
#-----------------------------------------------------------------------------------------------------------
# Lane-Emden Solutions
#-----------------------------------------------------------------------------------------------------------
n= 3.3 # Polytrope index
# dimensionless length : xi
# dimensionless density function : theta
# the gradient of theta : dtheta_dxi
# Empty arrays
xi=[]
theta=[]
dtheta_dxi=[]

# Delta x
dxi = 0.0001
# Boundary Conditions
theta.append(1.0)
dtheta_dxi.append(0.0) 
xi.append(10E-20) # Small non-zero value to avoid division by zero error

i=0 
while theta[i] > 0:  
    i=i+1
    dtheta_dxi_val = dtheta_dxi[i-1] - ( 2.0*dtheta_dxi[i-1]/xi[i-1] + (theta[i-1])**n)*dxi
    theta_val = theta[i-1] + dtheta_dxi_val*dxi
    xi_val = xi[i-1] + dxi
    # Append new calculated values to list. This way the last value of the list will be non-zero
    dtheta_dxi.append(dtheta_dxi_val)
    theta.append(theta_val)
    xi.append(xi_val)

# Stellar Surface values
xi_1 = xi[-1]
dtheta_dxi_1 = dtheta_dxi[-1]

# Calculate constants from Lane-Emden solutions
a_n = -(xi_1/3)*(1/dtheta_dxi_1)
b_n = 1/((n+1)*xi_1*(-1*dtheta_dxi_1))
c_n = (1/(4*np.pi*(n+1)))*(-1*dtheta_dxi_1)**2

Roh_c = a_n*3*Mstar/(4*np.pi*Rstar**3)
P_c = c_n*(G*Mstar**2)/(Rstar**4)
T_c = b_n*(G*Mstar*mu*m_H)/(Rstar*k_B)

K = ((4*np.pi)**(1/n)) * (G/(n+1)) * (xi_1**(-(n+1)/n)) * (-dtheta_dxi_1**((1-n)/n)) * (Mstar**((n-1)/n)) * (Rstar**((3-n)/n)) 
K = K.real
print(Roh_c, P_c, T_c, K)

#-----------------------------------------------------------------------------------------------------------
# Polytropic model
#-----------------------------------------------------------------------------------------------------------

alpha = np.sqrt(((n+1)*K*Roh_c**((1/n)-1))/(4*np.pi*G))
radius = []
radius[:] = [x * alpha for x in xi] # List of radius from xi
pressure = []
for i in theta:
    p = P_c*i**(n+1)
    pressure.append(p)

plt.plot(radius, pressure)
plt.show()
