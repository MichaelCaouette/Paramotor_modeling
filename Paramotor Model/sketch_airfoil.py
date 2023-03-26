# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 09:34:46 2022

Goal:
    Sketch the airfoil

@author: mcaoue2
"""


import numpy as np




class SketchAirFoil():
    def __init__(self):
        return
    
    
    def get_contour(self, angle=0, N_pts=100):
        
        # From https://www.researchgate.net/publication/312222678_Simple_analytic_equation_for_airfoil_shape_description
        B = 2
        T = 0.15
        C = 0.08
        P = 1
        E = 1
        R = 0
        
        # Parametric parameter
        t = np.linspace(0.001, 2*np.pi-0.001, N_pts)
        
        x = 0.5 + 0.5*np.abs(np.cos(t))**B/np.cos(t)
        
        y = (0.5*T*np.abs(np.sin(t))**B/np.sin(t)*(1-x**P) + 
             C*np.sin(x**E*np.pi) + 
             R*np.sin(x*2*np.pi) ) 
        
        # Rotate
        rot_M = np.matrix([[np.cos(-angle), -np.sin(-angle)], 
                           [np.sin(-angle),  np.cos(-angle)] ])
        
        cx, cy = np.dot(rot_M, [x, y])
        
        return np.array(cx)[0], np.array(cy)[0]
    
    
if __name__ == '__main__':  
    # Verify some property of the airfoil
        
    # =============================================================================
    #     Visualize the skecth
    # =============================================================================
    
    import matplotlib.pyplot as plt    


    self = SketchAirFoil()
    
    cont_x, cont_y = self.get_contour(angle=0.1)
    
    plt.figure()
    plt.plot(-cont_x, cont_y, color='r')
    plt.axis('equal')













