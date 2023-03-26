# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 14:34:33 2022

Goal:
    Force equation for the airfoil
Useful source:
    https://www.engineersedge.com/fluid_flow/aerodynamics_airfoil_theory_equations__15301.htm
    https://www.engineersedge.com/fluid_flow/long_symmetrical_airfoil_14046.htm
    
    Daytona wing
    https://www.itv-wings.com/fr/voiles-parapente-paramoteur/daytona.html
    
@author: mcaoue2
"""


import numpy as np


class ForceAirFoil():
    """
    Class for defining the airfoil and the force acting on it. 
    
    External paramter for determining the dynamic:
        Angle of attvk
        Velocity
    """
    
    def __init__(self, g=9.8067, air_density=1.225):
        """
        air_density:
            Density of the air (in kg/m^3)
            
        g:
            m/s^2
            Gravitational acceleration            
        """
        self.rho = air_density # kg/m^3
        self.acc_g = g
        
        # Set some physical quantity
        # From my Daytona wing, see the technical details here: (I have the 26)
        # https://www.itv-wings.com/fr/voiles-parapente-paramoteur/daytona.html
               
        # Projected area pf the airfoil, as seen from above. 
        # Used for both the lift anf the drag force     
        self.Ap  = 26 # m^2 (From my Daytona wing 26, see link above)
        
        # Mass of the airfoil. 
        # Used for the dynamical equation (f=ma) and gravity effect (mg)
        self.mass = 5.8 # kg (From my Daytona wing 26, see link above)
        
        # Frontal area
        self.envergure = 11.4 # m (From my Daytona wing 26, see link above)
        self.hauteur = 0.35 # m (Approximation width of my wing)
        self.A_frontal = self.envergure*self.hauteur
        
        # Aspect ratio (ProjectedArea/ corde**2 = envergue**2/ProjectedArea)
        self.aspect_ratio = 4.95 #Unitless, From Daytona spec instead of using the value
        

        
        # Drag coefficient at zero lift
        # Value from https://www.grc.nasa.gov/www/k-12/airplane/shaped.html#:~:text=The%20drag%20coefficient%20is%20a,of%20the%20velocity%20V%20squared.
        self.drag_coef_zero_lift = 0.45 
        
        # Initiate
        # Set arbitraty dynamic variable
        self.set_velocity(7, -0.2)
        self.set_AOA(10*np.pi/180)    
    
    def set_velocity(self, vx, vy):
        """
        In m/s. 
        Velocity of the airfoil. 
        (Such that the darg points in -vx, -vy)
        
        """
        self.vx, self.vy = vx, vy
        
        # Update the amplitude and angle of the velocity        
        self.v_square = self.vx*self.vx + self.vy*self.vy
        self.v_angle = np.angle(self.vx + 1j*self.vy)
        
    def set_AOA(self, angle_attack):
        """
        Set the angle of attack (AOA) in Radian.
        IT IS WITH RESPECT TO THE SURRUNDING LINEAR air FLOW 
        """
        self.angle_attack = angle_attack
        
    def _update_prms(self):
        """
        Update the parameters of the system.
        """
        
        # Velocity pressure
        self.VP = 0.5*self.rho*self.v_square   
 
        # Lift coefficient
        # Model from https://www.engineersedge.com/fluid_flow/aerodynamics_airfoil_theory_equations__15301.htm
        # There may exist a better model
        self.CL = self.get_lift_coef(self.angle_attack)

        
        # Drag coefficient
        # Model from https://www.engineersedge.com/fluid_flow/aerodynamics_airfoil_theory_equations__15301.htm
        # Induced drag from the lift (see https://www.grc.nasa.gov/www/k-12/airplane/dragco.html)
        # New URL https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-coefficient/
        self.c_induced  = self.CL*self.CL/(np.pi*self.aspect_ratio)
        # The total drag is the two
        self.Cd = self.drag_coef_zero_lift + self.c_induced
        
    def get_lift_coef(self, angle_attack):
        """
        Get the lift coefficient, which depends on the angle of attack
        """
        # Lift coefficient
        # Check here for some model VS angle of attack
        # http://www.aerospaceweb.org/question/aerodynamics/q0136.shtml
        # WHERE CAN I FIND A PROPER MODEL ? Try 0.3. Or solve it from the speed of my airfoil. 
        # self.k1 = 1.75         
        # # Angle de decrochage
        # # WHERE CAN I FIND A PROPER MODEL ? Try -20 deg
        # self.zero_lift_angle = 5*np.pi/180       
        # return self.k1*np.sin(angle_attack + self.zero_lift_angle)
        
        # Simple model to reproduce the figure of https://en.wikipedia.org/wiki/Lift_coefficient
        #Maximum value of the lift coefficient
        lift_max = 1.75
        # Angle for maximum lift
        angle_peak = 15*np.pi/180 
        # Negative angle for zero lift
        angle_zero = 5*np.pi/180
        # Stall angle maximum angle before stalling
        angle_stall = 25*np.pi/180
        # Put all this together
        cst = 0.5*np.pi/(angle_peak + angle_zero)
        curve = lift_max*np.sin((angle_attack+angle_zero)*cst)
        # For chopping the model
        u1 = np.heaviside(angle_attack+angle_zero, 0)
        u2 = np.heaviside(angle_stall- angle_attack, 0)
        
        return u1*u2*curve
       
    def lift_ampl(self):
        """
        Lift force amplitude 
        From the model on this page
        https://www.engineersedge.com/fluid_flow/aerodynamics_airfoil_theory_equations__15301.htm
        """
        
        # The lift
        return self.CL*self.VP*self.Ap
        # Direction is perp to v
        
    def drag_ampl(self):
        """
        Drag force amplitude.
        From the model on this page
        https://www.engineersedge.com/fluid_flow/aerodynamics_airfoil_theory_equations__15301.htm
        and
        https://www.engineersedge.com/fluid_flow/long_symmetrical_airfoil_14046.htm        
        """
        
        # The drag
        return self.Cd * self.A_frontal*self.VP
    
    def get_steady_state(self):
        """
        Get the steady state velocity of the airfoil. 
        This is when no external force is applied and the angle-of-attack is 
        constant. 
        
        I am using an analytical result that I derived from the toy model used 
        in this object (same model)

        Returns
        -------
        None.

        """
        # Update everything
        self._update_prms()       
        
        # A coefficient proportional to the lift force
        a = 0.5*self.rho*self.CL*self.Ap
        # A coefficient proportional to the drag force
        b = 0.5*self.rho*self.Cd*self.A_frontal
        
        v = (self.mass*self.acc_g/(a**2 + b**2)**0.5)**0.5
        
        angle = np.arctan(-b/a)
        
        return v, angle
        
                
        
    def get_all_force(self):
        """
        Get all force and direction. 
        Better, because it prevent to compute twice common parameters
        """
        # Update everything
        self._update_prms()
        
        # The drag force vector
        self.drag = self.drag_ampl()
        # It points opposite to the velocity of the airfoil 
        drag_x = -np.cos(self.v_angle)*self.drag
        drag_y = -np.sin(self.v_angle)*self.drag
        self.drag_vec = np.array([drag_x, drag_y])
        
        # The lift force vector
        self.lift = self.lift_ampl()
        # It is perpendicular to the velocity
        # Using these sines and cosines is as if we rotated the vector by 90 deg      
        # Proof: Make the dot product !
        lift_x = -np.sin(self.v_angle)*self.lift
        lift_y =  np.cos(self.v_angle)*self.lift
        self.lift_vec = np.array([lift_x, lift_y])

        # Gravitational force
        self.grav_vec = np.array([0, -self.mass*self.acc_g])

        # Summ
        return self.lift_vec, self.drag_vec, self.grav_vec
    
    def get_sum_force(self):
        """
        Get the some of all the force acting on the wing

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        list_forces = self.get_all_force()
        
        return np.sum(list_forces, axis=0)
        
        
        
        
        
if __name__ == '__main__':  
    # Verify some property of the airfoil
        
    # =============================================================================
    #     Check the lift coefficient VS angle of attack
    # =============================================================================
    
    import matplotlib.pyplot as plt
    
    list_angle = np.linspace(-15, 30, 100)*np.pi/180 #In rad
    
    self = ForceAirFoil()    
    list_cL =self.get_lift_coef(list_angle)
    
    
    plt.figure(tight_layout=True)
    plt.title('Lift coefficien')
    plt.plot(list_angle*180/np.pi, list_cL, label='My Airfoil',  linewidth=4)
    plt.xlabel('Angle of attack (deg)')
    plt.ylabel('Cl (Unitless)')
    plt.legend()        
    
    
    # =============================================================================
    #     Check the various forces
    # =============================================================================
            
    lift, drag, grav = self.get_all_force()
    F_total = self.get_sum_force()
    
    from sketch_airfoil import SketchAirFoil
    wing_sketch = SketchAirFoil().get_contour(angle=self.angle_attack)
    
    plt.figure(tight_layout=True)
    
    # Show the wing
    plt.plot(wing_sketch[0], wing_sketch[1], color='k', label='wing',  linewidth=4)    
    # Show the forces
    origin = [0.4, 0]
    scaling = 1/3e3
    plt.arrow(*origin, *lift*scaling, color='C1', label='Lift force')
    plt.arrow(*origin, *drag*scaling, color='C2', label='Drag force')
    plt.arrow(*origin, *grav*scaling, color='C3', label='Grav force')
    # plt.quiver(*origin, *lift, color='C1', label='Lift force', units='xy' ,scale=0.01)
    # plt.quiver(*origin, *drag, color='C2', label='Drag force', scale=5, scale_units='inches')
    # plt.quiver(*origin, *grav, color='C3', label='Grav force', scale=5, scale_units='inches')
    
    # Beautification
    plt.title('Lift coefficient')
    plt.legend()  
    plt.axis('equal')
    
    # =============================================================================
    # Check the various force as we change the Angle of Attack
    # =============================================================================

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
        
        
    # List all the angle to check up
    list_angle = np.linspace(-45, 45, 200)*np.pi/180
    # Calculate all the forces for each angle
    list_Flift = np.zeros(len(list_angle))
    list_Fdrag = np.zeros(len(list_angle))
    list_Fgrav = np.zeros(len(list_angle))
    list_Flift_angle = np.zeros(len(list_angle))
    list_Fdrag_angle = np.zeros(len(list_angle))
    list_Fgrav_angle = np.zeros(len(list_angle))
    for i, angle in enumerate(list_angle):
        self.set_AOA(angle) # Set the angle
        # Get the x,y coordinate of each forces
        Fl, Fd, Fg = self.get_all_force()  
        # Convert into magnitude and phase
        list_Flift[i], list_Flift_angle[i] = get_mag_angle(Fl[0], Fl[1])
        list_Fdrag[i], list_Fdrag_angle[i] = get_mag_angle(Fd[0], Fd[1])
        list_Fgrav[i], list_Fgrav_angle[i] = get_mag_angle(Fg[0], Fg[1])       
    
    # Plot all these results
    list_y = [list_Flift, list_Flift_angle*180/np.pi, 
              list_Fdrag, list_Fdrag_angle*180/np.pi, 
              list_Fgrav, list_Fgrav_angle*180/np.pi]
    list_label=['Lift (N)', 'Lift angle (deg)', 
                'Drag (N)', 'Drag angle (deg)',
                'Gravity (N)', 'Gravity angle (deg)']
    
    fig, axs = plt.subplots(nrows=len(list_y), ncols=1, 
                            sharex=True)#, tight_layout=True)
    for i, ax in enumerate(axs):
        ax.plot(list_angle*180/np.pi, list_y[i], color='C%d'%i)
        
        ax.legend([list_label[i]])
        ax.set_ylabel(list_label[i])
        
    ax.set_xlabel('Angle of attack (deg)')
    axs[0].set_title('Forces on the Airfoil')
    
    # =============================================================================
    #     Study the steady state (constant velocity)
    # =============================================================================
    
    # Here: Keep pursuing the analysis: Vary the angle of attack
    #       And set the speed to be this speed to compute the resulting
    #       lift and drag. Compare with 100kg*10m/s^2 = 1000N. 
    #       Is the lift force consistent with mg ? Does it have the right 
    #       amplitude to lift ?
    
    self.set_AOA(5*np.pi/180)
    self.mass = 100 # Kg, set the mass to something closer to my paramotor system
    v_steady, angle_steady = self.get_steady_state()
    print('v_steady (km/h), angle_steady (deg) = ', 
          v_steady*3.6, angle_steady*180/np.pi)
    print('This correspond, in m/s of vx, vy = ',
          v_steady*np.cos(angle_steady), 
          v_steady*np.sin(angle_steady))












