# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 15:10:18 2023

Goal:
    Create some GUI to monitor some variables during the simulation/game

@author: micha
"""

import matplotlib.pyplot as plt

import traceback
_p = traceback.print_last #Very usefull command to use for getting the last-not-printed error


class GUIMonitorVariables():
    """
    GUI to show a list of array. 
    It is meant to be uptade dynamically.
    """
    def __init__(self, x_axis, x_axis_label,
                   list_y_axis, list_label, title=''): 
        """
        Famous initiate function.

        x_axis:
            1D array (of float)
            x-axis shared with all the y-axis
        x_axis_label:
            string
            Label for the x_axis
        list_y_axis:
            list of 1D array (of float)
            Each element is an array to be plotted
        list_label:
            list of string
            Each element is the string used to label each element in 
            list_y_axis. 
        
        title='': (optional)
            String
            Text shown at the top of the GUI

        """    
        
        # # Initialize the figure and axis
        self.fig, self.axs = plt.subplots(nrows=len(list_y_axis), 
                                ncols=1, 
                                sharex=True,
                                tight_layout=True)
        
        # Update the plot
        self.update(x_axis, x_axis_label, list_y_axis, list_label,
                    title=title)
        
        self.fig.show()
        
        return 


    def update(self, x_axis, x_axis_label,
                   list_y_axis, list_label,title=''):
        """
        Update the plot
        
        x_axis:
            1D array (of float)
            x-axis shared with all the y-axis
        x_axis_label:
            string
            Label for the x_axis
        list_y_axis:
            list of 1D array (of float)
            Each element is an array to be plotted
        list_label:
            list of string
            Each element is the string used to label each element in 
            list_y_axis. 
            
        title='': (optional)
            String
            Text shown at the top of the GUI

        
        """
        
        # Put the array in each axes
        for i, ax in enumerate(self.axs):
            # Clear it first
            ax.cla()
            
            # Add the stuff
            ax.plot(x_axis, list_y_axis[i], color='C%d'%i)
            
            ax.legend([list_label[i]])
            
            ax.set_ylabel(list_label[i])
            
        ax.set_xlabel(x_axis_label)
        
        self.axs[0].set_title(title)  
        
        # Important. 
        #The following update the plot. 
        self.fig.canvas.draw_idle()  

if __name__=="__main__":
    # The the methods

    
    # Put some stuff in
    import numpy as np
    x = np.linspace(-1, 1, 300)
    ys = []
    str_ys = []
    for i in range(4):
        ys.append(i*np.sin(2*x*i))
        str_ys.append('Function %d'%(i+1))
 
    self = GUIMonitorVariables(x, 'Dummy x', ys, str_ys)
        
    # Do it again to test if it refresh properly with update
    import time
    print('Wait before refreshing...')
    time.sleep(3)
    print('Now refreshing')
    x = np.linspace(-1, 1, 300)
    ys = []
    str_ys = []
    for i in range(4):
        ys.append((4-i)**2*np.cos(2*x*i))
        str_ys.append('Function %d'%(i+1))
    
    self.update(x, 'Dummy x', ys, str_ys, title='Test')
        
    
    
        
    
    