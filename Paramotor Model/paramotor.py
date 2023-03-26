# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 11:56:35 2022

Goal:
    Descibre the force acting on the whole paramotor. Combining the airfoil
    and the motor. 

@author: mcaoue2
"""

import numpy as np

from airfoil import ForceAirFoil
from motor import ForceMotor
from solve_twomasses_reducedfreedom import TwoMassOneString


class ForceParamotor():
    """
    Class for defining the whole paramotor, which is a combination of the
    airfoil and the motor. 
    The combination is made such that they are "locked" together. 
    How ? By a tension from a connecting string. (the suspentes in real life)
    The model is a system of two masses (the wing and the motor) attached by 
    a string. Each of the masses feels a different total force PLUS the tension
    in the string. 
    
    I developped a generic solution for this type of system, for which I apply 
    to this particular system. 
       
    """
    
    def __init__(self, 
                 f_thrust,
                 g=9.8067, 
                 air_density=1.225, 
                 m_motor=100):
        """
        TODO add brake force VS time
        
        f_thrust:
            N
            Amplitude of the thrust force from the motor. 
            The signature is (t), where t is the time. (So it can vary with 
            time). And it must return a single number (the thrust force)
            
        
        air_density:
            Density of the air (in kg/m^3)
            
        g:
            m/s^2
            Gravitational acceleration            
        m:
            kg
            Mass of the motor + passager
            
        

        """
        self.f_thrust = f_thrust
        
        self.acc_g = g
        self.rho = air_density # kg/m^3
        self.m_motor = m_motor
        
        # Lenght of the line joining the wing to the motor
        # I would have to measure the lenght of my Daytona wing, because 
        # I don't find any official lenght here
        # https://www.itv-wings.com/fr/voiles-parapente-paramoteur/daytona.html
        self.dist = 2.5 # meter, approximate.
        
        # Define the components of the paramotor
        self.airfoil = ForceAirFoil(g=self.acc_g, 
                                    air_density=self.rho)
        self.motor = ForceMotor(g=self.acc_g, m=self.m_motor)
        
        # Equilibrium angle of attack of the wing. 
        # I need to check, but I am pretty sure that we don't fly at 0 angle
        # of attack
        self.AOA_offset = 10*np.pi/180 # In rad

        # =============================================================================
        # Set the physical system
        # =============================================================================
        # The physical system consist of two masses (the motor and the wing) 
        # connect by a string (the lines). 
        self.system = TwoMassOneString(m1=self.motor.mass,
                                       m2=self.airfoil.mass,
                                       distance=self.dist,
                                       F1 = self.F_motor,
                                       F2 = self.F_wing)
        
    def F_motor(self, t, vx, vy, phi, vphi):
        """
        Define the function for the force acting on the motor (masse 1),
        based on the dynamical variables of the system. 
        
        There are more input parameters than used. This is to match the 
        signature requested by the solver.
        """        
        
        # Set the thrust force
        self.motor.set_thrust(self.f_thrust(t))
        # Set the angle of the force
        angle_thrust = np.pi/2 - phi
        self.motor.set_angle_thrust(angle_thrust)
        
        return self.motor.get_sum_force()
        
        
    def F_wing(self, t, vx, vy, phi, vphi):
        """
        Define the function for the force acting on the aifoil (masse 2),
        based on the dynamical variables of the system. 
        
        There might be more input parameters than used. This is to match the 
        signature requested by the solver.
        
        vx,vy is the speed of what ?T The airfoil or the motor ?
        The motor I think !
        """
        
        # Set the parameters that defines the forces
        # Angle of attack # \m/
        # It depends on the relative position with the motor AND it is with 
        # respect to the flow. 
        # First the angle of the airfoil
        beta = np.pi/2 - phi
        bx, by = np.cos(beta), np.sin(beta) # Unit vector 
        # The angle of attack is the angle of the airfoil with respect to the 
        # wind flow (in the negative direction). Or equivalently, the airfoil 
        # velocity vector
        # The airfoil velocity is the sum of the bottom speed and the angular
        # velocoty
        vx_airfoil = -self.dist*np.sin(phi)*vphi + vx
        vy_airfoil = +self.dist*np.cos(phi)*vphi + vy
        # I am using a general formula, because I want the full range between
        # -pi to +pi. Just using the sin or cos definition restricts the range
        
        # We find the cos and sin of the angle between the two vector
        # MODULO the magnitude, which is not important to get the angle
        cosy = (bx*vx_airfoil + by*vy_airfoil) # From vector dot   product rule
        siny = (bx*vy_airfoil - by*vx_airfoil) # From vector cross product rule
        # Retrieve the angle by using the complex plane methode (easy)      
        AOA = np.angle(cosy + 1j*siny) + self.AOA_offset
        self.airfoil.set_AOA(AOA)       
        # The velocity of the airfoil
        self.airfoil.set_velocity(vx_airfoil, vy_airfoil)
        
        return self.airfoil.get_sum_force()
    
    def update_system(self, t, vx, vy, phi, vphi):
        """
        Update all the internal variable of the subssystems. 
        Return nothing
        
        """
        # Update the wing
        self.F_wing(t, vx, vy, phi, vphi)
        # Update the motor
        self.F_motor(t, vx, vy, phi, vphi)
        # All done
        return
        
        
    
    def solve_dynamic(self, list_t, init_parms):
        """
        EXPLANATION PLEASE
        
        init_parms:
            The initial condition of the 
            degree of freedoms of the system. 
            AKA: The velocity vector of the motor
            and the angle and angular velocity of the line. 
            ['vx1','vy1', 'phi', 'vphi']
        """
        

        # sovle it
        self.list_var, self.sol = self.system.solve_dynamic(list_t, 
                                                            init_parms,
                                                            method='RK45',
                                                            rtol=1e-4,
                                                            atol=1e-8)
        
        return self.list_var, self.sol

    def solve_iteratively(self,
                          init_parms, 
                          dt, t0, duration):
        """
        EXPLAIN MORE PLEASE
        
        Solve the system iteratively. 

        init_parms:
            The initial condition of the 
            degree of freedoms of the system. 
            AKA: The velocity vector of the motor
            and the angle and angular velocity of the line. 
            ['vx1','vy1', 'phi', 'vphi']     
        dt:
            (float)
            Time step for the solver
        t0:
            (float)
            Initial time.
        duration:
            (float)
            Total time for solving the system
            
        
        Returns
        -------
        None.

        """
        # My goal is to make the code general for any number of variable. 

        # For this specific system it should be [vx1, vy1, phi, vphi]
        list_var = init_parms # Initiate the variables
        N_var = len(list_var) # Number of variable to be solved
        
        # Initiate the array of variables that we save
        self.list_t = []
        self.list_vars = [] # List of list of the variablesé 
        for ii in range(N_var):
            # Each element of self.list_vars will be an array for the variable
            # number ii
            self.list_vars.append([])   
        
        # We will now solve the system at each increment of time. 
        t_acc = t0
        t_end = t0 + duration
        while (t_acc < t_end):
            # Obtain the time derivative
            out = self.system.derivatives(t_acc, list_var)
            list_derivative = out            
            # Update each variable
            for ii in range(N_var):
                # Solve with a simple Newton method for now
                list_var[ii] += list_derivative[ii]*dt
                # Save the value
                self.list_vars[ii].append(list_var[ii])
                
            # Save the time 
            self.list_t.append(t_acc)
            
            # Increment the total time elapsed for the next iteration
            t_acc += dt
        # convert into numpy array (so the user can do math !)
        return np.array(self.list_t), np.array(self.list_vars)
        
    def solve_quick(self,
                          init_parms, 
                          dt, t0, duration):
        """
        Function to quickly iterate and outputting only the last variable. 
        Usefull for game dynamic.
        
        Solve the system iteratively. 

        init_parms:
            The initial condition of the 
            degree of freedoms of the system. 
            AKA: The velocity vector of the motor
            and the angle and angular velocity of the line. 
            ['vx1','vy1', 'phi', 'vphi']     
        dt:
            (float)
            Time step for the solver
        t0:
            (float)
            Initial time.
        duration:
            (float)
            Total time for solving the system
            
        
        Returns
        -------
        None.

        """
        # My goal is to make the code general for any number of variable. 

        # For this specific system it should be [vx1, vy1, phi, vphi]
        list_var = init_parms # Initiate the variables
        N_var = len(list_var) # Number of variable to be solved
        
        
        # We will now solve the system at each increment of time. 
        t_acc = t0
        t_end = t0 + duration
        while (t_acc < t_end):
            # Obtain the time derivative
            out = self.system.derivatives(t_acc, list_var)
            list_derivative = out            
            # Update each variable
            for ii in range(N_var):
                # Solve with a simple Newton method for now
                list_var[ii] += list_derivative[ii]*dt
            
            # Increment the total time elapsed for the next iteration
            t_acc += dt
        # convert into numpy array (so the user can do math !)
        return list_var    
        
        


if __name__ == "__main__":
    """
    The following code is executed only when this script is run alone. 
    Used to test the code
    """        
    import matplotlib.pyplot as plt 
    
    def my_thrust(t):
        # constant force
        return 0
    
    # Initial condition
    # Its in the order ['vx1', 'vy1', 'phi', 'vphi']
    # In the air, with some existing velocity
    init_var = [ 7,    # From the steady state analysis
                -0.47, # From the steady state analysis
                1.0*np.pi/2+0.0,
                0.0]
    # Time axis, in sec
    N_pts = 300
    t = np.linspace(0, 1, N_pts)*0.9
    
    # TO DO, various initial angle for taking off and compare how it helps
    
    self = ForceParamotor(my_thrust)
    
    print('Solving the dynamic...')
    z, ss = self.solve_dynamic(t, init_var)
    print('Done')
    print('Solve again with simple iterative method...')
    t_iter, z_iter = self.solve_iteratively(init_var, 
                                            dt=0.001, 
                                            t0=t[0], 
                                            duration = t[-1]-t[0])
    print('Done')
    
    
    list_label=['vx1 (m/s)', 
                'vy1 (m/s)',
                'phi (rad)', 'vphi (rad/s)']
    # Plot that up
    fig, axs = plt.subplots(nrows=len(init_var), 
                            ncols=1, 
                            sharex=True)#, tight_layout=True)
    for i, ax in enumerate(axs):
        # Solution of the solver
        ax.plot(t, z[i], color='C%d'%i)
        # Solution of the simple iterative method
        ax.plot(t_iter, z_iter[i], 'k--')
        ax.legend([list_label[i],'Simple iterative method'])
        
        ax.set_ylabel(list_label[i])
        
    ax.set_xlabel('t (s)')
    axs[0].set_title('System')
    
    # Check how other variable evolved during the simulation
    list_AOA = [] # Angle of attack
    list_lift_coef = [] # Lift coefficient
    for i in range(N_pts):
        # Set the force, this will update the variables in the solver
        vx, vy, phi, vphi = z.T[i] # Get the variable at this iteration
        self.update_system(t[i], vx, vy, phi, vphi)
        # Update the variable that we care about
        aa = self.airfoil.angle_attack
        list_AOA.append(aa)
        list_lift_coef.append(self.airfoil.get_lift_coef(aa))
    # Make that into numpy array
    list_AOA = np.array(list_AOA)
    list_lift_coef = np.array(list_lift_coef)

    # Plot that up
    # Note to myself: I should make this a function because I use it way too 
    # much. 
    x_axis = t
    list_y_axis = [list_AOA*180/3.14159,
                   list_lift_coef]
    list_label  = ['Angle of Attack (deg)',
                   'Lift coefficient (Unitless)']
    
    fig, axs = plt.subplots(nrows=len(list_y_axis), 
                            ncols=1, 
                            sharex=True)#, tight_layout=True)
    for i, ax in enumerate(axs):
        ax.plot(x_axis, list_y_axis[i], color='C%d'%i)
        
        ax.legend([list_label[i]])
        
        ax.set_ylabel(list_label[i])
        
    ax.set_xlabel('t (s)')
    axs[0].set_title('System variables')    
 
    
    
        
        
        
        
        
        
        