# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 09:13:57 2022

Goal:
    Describe the force acting on the motor part of the paramotor

@author: mcaoue2
"""

import numpy as np


class ForceMotor():
    """
    Class for defining the motor and body mass and the force acting on it. 
    
    External paramter for determining the dynamic:
        Thrust force (From the back helix)
        Angle_thrust (Angle of the where the thrust force is applied)
    """
    
    def __init__(self, g=9.8067, m=100):
        """
        
        g:
            m/s^2
            Gravitational acceleration
            
        m:
            kg
            Mass of the motor + passager
        """
        # Physical property
        
        # Gravitational acceleration
        self.acc_g = g        
        # Mass of passager + motor
        self.mass = m # kg
        
        
    def set_thrust(self, thrust):
        self.thrust_force = thrust        
        
    def set_angle_thrust(self, angle_thrust):
        self.angle_thrust = angle_thrust
        
    def get_sum_force(self):
        """
        Get all the force acting on the system
        """
        
        # Trust force
        thrust_x = np.cos(self.angle_thrust)
        thrust_y = np.sin(self.angle_thrust)
        self.thrust_vec = np.array([thrust_x, thrust_y])

        # Gravitational force
        self.grav_vec = np.array([0, -self.mass*self.acc_g])
        
        # Sum
        return self.thrust_vec + self.grav_vec
        
        
        
if __name__ == '__main__':  
    # Verify some property
        
    # =============================================================================
    #     Visualize the forces
    # =============================================================================
    
    
    list_angle = np.linspace(-15, 30, 100)*np.pi/180 #In rad
    
    self = ForceMotor()    
    self.set_thrust(200)
    self.set_angle_thrust(np.pi/180*15)
    F_total = self.get_sum_force()
    print('F_total = ', F_total)
    

        


        
        
        