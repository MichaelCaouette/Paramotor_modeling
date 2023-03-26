# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 09:29:23 2022

goal:
    Solve the dynamic of the physical model: two masses attached by a string. 

@author: client
"""


from scipy.integrate import solve_ivp
import numpy as np

class TwoMassOneString():
    """
    Solver for a specific physical system. 
    For describnig the system and solving its dynamic.
    
    THIS SOLVER IS SIMPLIFIED FOR THE PARAMOTOR AIRCRAFT. It just take, as degree of 
    freedom, the parameters in the equation of motion present in the aircraft.
    (Velocity of one point and angle)
    
    Physical model for two point masses attached by a string of fixed lenght. 
    Each masses feels an independent forces, plus the tension in the string. 
    The lenght of the string is fixed and the string is straight. Such that the
    distance between the two masses is fixed. 
    This extract constrain reduces the number of degree of freedom. 
    The degree of freedom that we solve are
        - vx, vy --> The velocity of the first mass
        - phi    --> The angle of the string
        - vphi   --> Time derivative of the string 
                        (impacts the tension in the string and velocity of m2)
    
    Note that we are not solving for the spatial position of the masses. This 
    is because we assum that the external force does not depend on position, 
    which is the case of the aircraft (for which this class is created).
    The solver that consider all these freedom, use the script 
    'solve_twomasses_allfreedom.py'. 
    
    """
    def __init__(self, 
                 m1, m2, distance,
                 F1, F2):
        """
        m1, m2:
            Kg
            Masses of the first and second point masses
        distance:
            m
            Distance (fixed) between the two masses. 
        F1, F2:
            N
            Function describing the forces acting on masse 1 and masse 2. 
            The signature of these function must be 
            (t, vx, vy, phi) and it must return (Fx, Fy) 
            (ie the two components)
            Note that:
                phi is the angle of the vector
                r2-r1 (AKA the string joining the two masses)
            
        """
        self.m1, self.m2 = m1, m2
        self.dist = distance
        self.F1, self.F2 = F1, F2 
        
    def derivatives(self, t, list_var):
        """
        Compute the set of 1rst order diff equation describing the dynamic.
        
        We use the 4 degree of freedom related to the first mase and the angle 
        of the string:
            list_var = vx1, vy1, phi, vphi
            Where:
                vx1, vy1 is the velocity of this first mass
                phi and vphi is the angle and time derivative of the angle of 
                the string (angle of the vector r2-r1)
        
        Return:
            list of the 4 time derivatives
        """
        # We have 4 dynamical varibales:
        vx1, vy1, phi, vphi = list_var
        
        # Some parameters to simplify the equations latter on
        cosy, siny = np.cos(phi), np.sin(phi)
        vx2 = -self.dist*siny*vphi + vx1
        vy2 =  self.dist*cosy*vphi + vy1
        
        # The forces acting on each masses
        force_1 = self.F1(t, vx1, vy1, phi, vphi) # It is a vector (fx, fy)
        force_2 = self.F2(t, vx2, vy2, phi, vphi) # It is a vector (fx, fy)
        
        # Simplification
        qq = 1/(self.m2*self.dist)
        ratio = self.m2/self.m1
        beta = 1/(1 + ratio)
        a = force_2[0] - ratio*force_1[0]
        b = force_2[1] - ratio*force_1[1]
        
        # Compute the tension
        # Amplitue
        T = beta*(b*siny + a*cosy + self.m2*self.dist*vphi*vphi)
        # Vector of the tension force
        tension = [T*cosy, T*siny]
        
        # Associate the time derivatives with the proper variables
        # The velovities
        d_phi = vphi
        # The second derivatives (Derivative of velocities)
        d_vphi = qq*(b*cosy - a*siny)
        d_vx1 = (force_1[0] + tension[0])
        d_vy1 = (force_1[1] + tension[1])   
        
        return (d_vx1, d_vy1, d_phi, d_vphi)
            
    def solve_dynamic(self, list_t, y0, **kwargs):
        """
        EXPLAIN
        """
        self.sol = solve_ivp(self.derivatives,
                             [min(list_t), max(list_t)], 
                             y0,
                             dense_output=True,
                             **kwargs)        
        list_var = self.sol.sol(list_t)
        return list_var, self.sol
    

if __name__ == "__main__":
    """
    The following code is executed only when this script is run alone. 
    """    
    import matplotlib.pyplot as plt 
    
    # Test with known solution system. 
    # =============================================================================
    # A simple pendulum !    
    # =============================================================================
    # To keep m2 fixed in space, the force F2 must be equally the opposite of
    # the tension on the string. However, we don`t know this tension a priori. 
    # Therefore I will, instead, make F2 to be an upward force and m2 >> m1
    m1, m2 = 3, 9875
    g = 9.81
    dist = 0.5
    # Initial condition
    # Its in the order ['vx1', 'vy1', 'phi', 'vphi']
    init_var = [0, 0,
                np.pi/2-0.2, 
                0.0]
    
    def F1(t,  vx, vy, phi, vphi):
        # The force is actually constant and doesn`t depend on the variable. 
        return [0, -m1*g]
    
    def F2(t, vx, vy, phi, vphi):
        # We simulate a pull in the upward direction
        # The pull is equal to the total weight of the system minus the weight
        # of m2. Therefor just m1*g
        return [0, m1*g]        
    
    # Time axis
    t = np.linspace(0, 5, 300)
    
    # Expected solution
    exp_w = (g/dist)**0.5
    A = np.pi/2 - init_var[2]
    B = (-init_var[3])/exp_w
    exp_theta = A*np.cos(exp_w*t) + B*np.sin(exp_w*t)
    exp_phi = np.pi/2 - exp_theta
    
    # Solve the system !
    self = TwoMassOneString(m1, m2, 
                            dist,
                            F1, F2)
    
    z, ss = self.solve_dynamic(t, init_var)
    
    list_label=['vx1 (m/s)', 
                'vy1 (m/s)',
                'phi (rad)', 'vphi (rad/s)']
    # Plot that up
    fig, axs = plt.subplots(nrows=len(init_var), ncols=1, 
                            sharex=True, tight_layout=True)
    for i, ax in enumerate(axs):
        ax.plot(t, z[i], color='C%d'%i)
        if i==2:
            ax.plot(t, exp_phi, '--k')
            ax.legend([list_label[i], 'Harmonic Oscillator'])
        else:
            ax.legend([list_label[i]])
        
        ax.set_ylabel(list_label[i])
        
    ax.set_xlabel('t (s)')
    axs[0].set_title('System')
 
    
    #TODO Plot phi and expected PHI
    
    
    # =============================================================================
    # Driven pendulum !  
    # =============================================================================
    # I think that we can have an mostly analytic solution to this one. 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

