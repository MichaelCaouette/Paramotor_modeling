# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 17:55:55 2023

Goal:
    Test the aifoil object.
    Test some prediction of the model. 

@author: micha
"""
import numpy as np
import matplotlib.pyplot as plt

from airfoil import ForceAirFoil

def get_mag_angle(fx, fy):
    """
    Get the magnitude and angle of a vector (fx, fy)

    Parameters
    ----------
    fx : TYPE
        DESCRIPTION.
    fy : TYPE
        DESCRIPTION.

    Returns
    -------
    mag, angle
    The magnitude and angle of the input vector. The angle is in radian

    """
    mag = (fx*fx + fy*fy)**0.5
    # Retrieve the angle by using the complex plane methode (easy)      
    angle = np.angle(fx + 1j*fy) 
    
    return mag, angle
    
# =============================================================================
#     Study the steady state (constant velocity)
# =============================================================================

my_mass = 100 # kg, mass applied to the wing

self = ForceAirFoil()

# Cheat a bit for setting the mass (should be from a method !)
self.mass = my_mass

# Quick test/check up to see if the value makes sense
self.set_AOA(2*np.pi/180)
v_steady, angle_steady = self.get_steady_state()
# For the regime of my wing and weight, I expect a speed of 
# few tens of km/h and an angle between -90 deg and 0 deg. 
print('v_steady (km/h), angle_steady (deg) = \n', 
      v_steady*3.6, angle_steady*180/np.pi)

# We will now vary the angle of attack and check the steady state speed
# Let's not exceed the stall angles
list_AoA = np.linspace(-4.9, 24.9, 500)*np.pi/180
# Prepare the list to study
list_v, list_av = np.zeros((2, len(list_AoA)))
list_Flift = np.zeros(len(list_AoA))
list_Fdrag = np.zeros(len(list_AoA))
list_Fgrav = np.zeros(len(list_AoA))

for i, angle in enumerate(list_AoA):
    self.set_AOA(angle) # Set the angle
    list_v[i], list_av[i]= self.get_steady_state()
    # Set the steady state velocity and compute the the forces
    self.set_velocity(list_v[i]*np.cos(list_av[i]),
                      list_v[i]*np.sin(list_av[i]) )
    # Get the x,y coordinate of each forces
    Fl, Fd, Fg = self.get_all_force()  
    # Convert into magnitude and phase
    list_Flift[i], _ = get_mag_angle(Fl[0], Fl[1])
    list_Fdrag[i], _ = get_mag_angle(Fd[0], Fd[1])
    list_Fgrav[i], _ = get_mag_angle(Fg[0], Fg[1])     

 
# Get the cartesian speed
list_vx =  list_v*np.cos(list_av)  
list_vy =  list_v*np.sin(list_av)   
# Get the AoA for which the speed is maximal
i_max = np.argmax(list_v) 
AoA_max_v = list_AoA[i_max]
# And when the horizontal speed is maximal
i_max_x = np.argmax(list_vx)
AoA_max_vx = list_AoA[i_max_x]
# And when the vertical speed is maximal
i_max_y = np.argmax(list_vy)
AoA_max_vy = list_AoA[i_max_y]

# Test that the sum of force is actually zeros !
i_test = 5
self.set_AOA(list_AoA[i_test])
self.set_velocity(list_vx[i_test], list_vy[i_test])

print()
print('Test Sum of force: ', self.get_sum_force())
print()

# Plot the information !
plt.figure(tight_layout=True)
plt.subplot(311)
plt.plot(list_vx*3.6, list_vy*3.6)
plt.plot(list_vx[i_max]*3.6, list_vy[i_max]*3.6, 'or', 
         label='Maximum speed, AoA = %.1f deg'%(AoA_max_v*180/np.pi ))
plt.plot(list_vx[i_max_x]*3.6, list_vy[i_max_x]*3.6, 'ob', 
         label='Max horizontal speed, AoA = %.1f deg'%(AoA_max_vx*180/np.pi ))
plt.plot(list_vx[i_max_y]*3.6, list_vy[i_max_y]*3.6, 'og', 
         label='Max vertical speed, AoA = %.1f deg'%(AoA_max_vy*180/np.pi ))

plt.legend()
plt.xlabel('Vx (km/h)')
plt.ylabel('Vy (km/h)')
plt.title('Velocity as the angle of attack changes\n'+
          'Mass on the wing = %.1f kg'%my_mass)

# Show the speeds VS angle of attacks
plt.subplot(312)
plt.plot(list_AoA*180/np.pi, list_vy*3.6, label='Vy') 
plt.plot(list_AoA*180/np.pi, list_vx*3.6, label='Vx') 
# Show the cool points
plt.plot(list_AoA[i_max]*180/np.pi, list_vy[i_max]*3.6, 'or')
plt.plot(list_AoA[i_max]*180/np.pi, list_vx[i_max]*3.6, 'or')
plt.plot(list_AoA[i_max_x]*180/np.pi, list_vx[i_max_x]*3.6, 'ob')
plt.plot(list_AoA[i_max_y]*180/np.pi, list_vy[i_max_y]*3.6, 'og')
plt.legend()
plt.xlabel('Angle of attack (deg)')
plt.ylabel('Speeds (km/h)')
plt.title('Velocity VS angle of attack changes')


# Show the forces VS Angle of Attack
plt.subplot(313)
list_y = [list_Flift,  list_Fdrag, list_Fgrav]
list_label=['Lift',
            'Drag',
            'Gravity']

for i in range(len(list_y)):
    plt.plot(list_AoA*180/np.pi, list_y[i], 
             color='C%d'%i, label=list_label[i])
    # Also show the points
    plt.plot(list_AoA[i_max  ]*180/np.pi, list_y[i][i_max  ], 'or')
    plt.plot(list_AoA[i_max_x]*180/np.pi, list_y[i][i_max_x], 'ob')
    plt.plot(list_AoA[i_max_y]*180/np.pi, list_y[i][i_max_y], 'og')
             
plt.legend()
plt.ylabel('Amplitude of force (N)')
plt.xlabel('Angle of attack (deg)')
plt.title ('Forces on the Airfoil')








