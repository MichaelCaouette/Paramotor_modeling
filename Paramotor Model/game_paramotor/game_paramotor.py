# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 17:37:38 2023

Goal:
    Implement the paramotor equation to see how the parameters of the 
    aircraft evolve with time. 

@author: mcaoue2
"""

from paramotor import ForceParamotor

# Just type _p() in the console 
import traceback
_p = traceback.print_last

# Import the pygame module
import pygame

# Import random for random numbers
import random

# The legendary numpy
import numpy as np

# This is to show the evolution of some variables
from gui_monitor_variables import GUIMonitorVariables
# import matplotlib.pyplot as plt

# Import pygame.locals for easier access to key coordinates
# Updated to conform to flake8 and black standards
# By importing specific constants from pygame.locals, 
# you can use the syntax <CONSTANT> instead. 
# This will save you some keystrokes and improve overall readability.
from pygame.locals import (
    RLEACCEL,    
    K_UP,
    K_DOWN,
    K_LEFT,
    K_RIGHT,
    K_SPACE,
    K_ESCAPE,
    KEYDOWN,
    QUIT,
)

# =============================================================================
# Define constants
# =============================================================================
# Screen width and height
SCREEN_WIDTH  = 800
SCREEN_HEIGHT = 600
# The number of tiles in each direction
N_TILES_X = 10
N_TILES_Y = 9
# Variables used to tag various elements
TAG_CLOUD = 1
TAG_NOTHING = 0
# Frame per sec
FPS = 30 
# Conversion
deg2rad = np.pi/180
rad2deg = 180/np.pi

path_img_player = "paramotor_v1_shrinked.png"
path_img_cloud  = "cloud.png"

# =============================================================================
# Initiate Now. Pygame is useful for operations used in other initiation
# =============================================================================
# Initialize pygame
pygame.init()
# Create the screen object
# The size is determined by the constant SCREEN_WIDTH and SCREEN_HEIGHT
screen = pygame.display.set_mode((SCREEN_WIDTH, SCREEN_HEIGHT))


# =============================================================================
# Load the images for the tiles
# =============================================================================
cloud_img= pygame.image.load(path_img_cloud).convert()
cloud_img.set_colorkey((0, 0, 0), RLEACCEL)   
TILE_IMG_WIDTH, TILE_IMG_HEIGHT = cloud_img.get_size()


# =============================================================================
# Other constants that depend on the setting
# =============================================================================
# The size of the tiles (# of pixel)
TILE_WIDTH  = TILE_IMG_WIDTH  # Redondancy for now
TILE_HEIGHT = TILE_IMG_HEIGHT # Redondancy for now
# Map width and height (Maybe bigger than the view screen, or smaller)
MAP_WIDTH  = TILE_WIDTH  * N_TILES_X
MAP_HEIGHT = TILE_HEIGHT * N_TILES_Y
# How many tile to add such that they don't dispear before leaving the screen 
N_TILE_ADD_X = 1
N_TILE_ADD_Y = 1
# How many tiles fit in the screen in each direction
# Add some tile such that they don't disapear
N_TILE_TO_SHOW_X = int(round(SCREEN_WIDTH /TILE_WIDTH  ) ) + N_TILE_ADD_X*2
N_TILE_TO_SHOW_Y = int(round(SCREEN_HEIGHT/TILE_HEIGHT ) ) + N_TILE_ADD_Y*2
# Useful list
list_map_ix = list(range(N_TILES_X) )
list_map_iy = list(range(N_TILES_Y) )
# Time between frame 
time_btw_frame = 1000/FPS # millisecond
time_btw_frame_sec = time_btw_frame/1000 # In sec

# =============================================================================
# Generate the map
# =============================================================================
map_cloud = np.zeros((N_TILES_X, N_TILES_Y)) + TAG_NOTHING
prob_cloud = 5 # The likelihood to get a cloud on a tile is 1/prob_cloud
for i in range(N_TILES_X):
    for j in range(N_TILES_Y):
        # One chance over prob_cloud to have a cloud
        if random.randint(1, prob_cloud) == prob_cloud:
            map_cloud[i][j] = TAG_CLOUD  

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



# Define a Player object by extending pygame.sprite.Sprite
# The surface drawn on the screen is now an attribute of 'player'
class Player(pygame.sprite.Sprite):
    """
    The player as a velocity
    """
    def __init__(self):
        super(Player, self).__init__() # uses .super() to call the .__init__() method of Sprite
        # Instead of a surface, use an image for a better-looking sprite
        self.surf_ori = pygame.image.load(path_img_player).convert() # the .convert() call optimizes the Surface, making future .blit() calls faster.
        # The RLEACCEL constant is an optional parameter that helps pygame render more quickly on non-accelerated displays.
        # uses .set_colorkey() to indicate the color pygame will render as transparent. 
        self.surf_ori.set_colorkey((255, 255, 255), RLEACCEL)  
    
        # Center the player
        self.rect_ori = self.surf_ori.get_rect()
        self.rect_ori.centerx = (SCREEN_WIDTH -self.surf_ori.get_width() )/2
        self.rect_ori.centery = (SCREEN_HEIGHT-self.surf_ori.get_height())/2
        
        # Initiate the physical property
        # The thrust is set to match the Atom80 thrust force hereL
        # https://www.americanparagliding.com/impuls/atom80.htm
        self.F_trust_amplitude = 55*9.8 # Newton,  55kg * 9.8m/s**2 
        self.physical_system = ForceParamotor(f_thrust=self.model_thrust_force)
        # Initial dynamical variable
        # Its in the order ['vx1', 'vy1', 'phi', 'vphi']
        # In m/s for the speed
        self.list_dyn_var = [ 7,    # From the steady state analysis
                             -0.47, # From the steady state analysis
                             1.0*np.pi/2+0.0,
                             0.0]        
        
        # Indication of the state
        self.state_thrust = 0 # Either 0 or 1, for no thrust or full thrust
        self.state_brake = 0 # Either 0 or 1, for no brake or full brake. 
        
        # update everything (also initiate surf and rect)
        self._update_variables()

    
    def model_thrust_force(self, t):
        """
        Used to set the thrust force in the dynamic solver. 

        Parameters
        ----------
        t : float
            Dummy variable (not used) required to matches the 

        Returns
        -------
        Amplitude of the thrust force

        """
        
        return self.state_thrust*self.F_trust_amplitude
        
        
    def update(self, pressed_keys):
        # The pressed keys will update some of the physical property
        
        # The space sets the thrust force of the motor.
        # For now, because I don't see how to handle continuous variables,
        # it will be full throttle or no throttle at all. 
        if pressed_keys[K_SPACE]:
            self.state_thrust = 1
        else:
            # We stopped the thrust force
            self.state_thrust = 0 
            
        # The left key sets the brake
        if pressed_keys[K_LEFT]:
            self.state_brake= 1
        else:
            # We stopped the brake
            self.state_brake = 0 
        
        # To be added later:
        # - A key to set the trim 
        
        # update everything
        self._update_variables()
        
    def _update_variables(self):
 
        # Note to be deleted:
            # We might re-initiate the physic solver since we update the 
            # function for the thrust force
            
        # Update the physics
        # Use a quick-iterative-not-perfect method for now
        out = self.physical_system.solve_iteratively(self.list_dyn_var, 
                                                dt=time_btw_frame_sec/100, 
                                                t0=0, # Origine of time for the init variable
                                                duration = time_btw_frame_sec)
        # For now, we are not using the time form the output. 
        # Maybe at some point
        times, list_all_vars = out
        
        # Extract other quantity from the dynamic
        # It is a list of variable for each time step
        vx1_s, vy1_s, phi_s, vphi_s = list_all_vars
        # Just take the last element, which are the most up to date
        self.list_dyn_var = np.array([vx1_s[0], vy1_s[0], phi_s[0], vphi_s[0]])
        vx1, vy1, phi, vphi  = self.list_dyn_var
        
        # Debug
        print(vx1, vy1, phi, vphi)
        self.physical_system.update_system(0, vx1, vy1, phi, vphi)
        # Get the center of mass velocity
        # Speed and mass of the wing
        vx2 = self.physical_system.airfoil.vx
        vy2 = self.physical_system.airfoil.vy
        m2  = self.physical_system.airfoil.mass
        # For the motor, the speed is vx1 and vx2
        m1 = self.physical_system.motor.mass
        # Formula of CM velocity
        M = m1 + m2
        self.vx = (m1*vx1 + m2*vx2)/M
        self.vy = (m1*vy1 + m2*vy2)/M
        
        # The angle
        self.angle_system = np.pi/2 - phi
        
        # Update the image
        # Need to 
        # Set the surface and the rect to be at a certain angle
        self.surf = pygame.transform.rotate(self.surf_ori, 
                                            -self.angle_system *rad2deg)
        self.rect = self.surf.get_rect(center = self.rect_ori.center)
        
    def get_info(self):
        """
        Return a lot of local variable about the system. 
        Useful for keeping informed about the state of the system

        Returns
        -------
        self.list_label_info, self.list_info

        """
        
        vx1, vy1, phi, vphi  = self.list_dyn_var        
        self.list_info = [vx1,
                          vy1,
                          phi*rad2deg,
                          vphi*rad2deg,
                          self.vx,
                          self.vy,
                          self.model_thrust_force(0),
                          self.state_brake]
        # A string for each stuff saved
        self.list_label_info = ['vx1',
                                'vy1',
                                'phi (deg)',
                                'vphi (deg/sec)',
                                'vx',
                                'vy',
                              'thrust_force',
                              'state_brake\n(not used yet)']
        
        return self.list_label_info, self.list_info
        
        
                


def circularList(list_original, i_start, N_list):
    # Create a list of element from "list_original" that circulate periodically
    # from the "i_start" and N_list time. 
    # It is used to generate a periodic boundary condition on the map
    # More explanation later, but it works !
    
    # Lenght of original list tells us how to modulate
    N_modulation = len(list_original)
    
    # The initial index in the original list as to be modulated
    imin = i_start%N_modulation
    my_list = []
    for j in range(N_list):
        # The new index is shifted and modulated
        i_new = (imin+j)%N_modulation
        my_list.append(list_original[i_new])
    return my_list
        
            
# Setup the clock for a decent framerate
clock = pygame.time.Clock()

# Instantiate player. 
player = Player()

# Create groups to hold cloud sprites, and all sprites
# - clouds is used for position updates
# - all_sprites is used for rendering
# clouds = pygame.sprite.Group()
all_sprites = pygame.sprite.Group()
all_sprites.add(player) # Add the player

# Variable to keep the main loop running
running = True

# The following variable will need to be put into a class later
# The player starts at the center of the map
# It is in the map coord system
player_x, player_y = MAP_WIDTH/2, MAP_HEIGHT/2


# =============================================================================
# Variable to be saved
# =============================================================================
want_save_data = True
# Maximum element to save for each variable It will overwrite the first element
N_max_save = 500
# List that, which will simplify how we loop over 
save_label, saved_infos = player.get_info()
list_saved_things = []
for i in range(len(saved_infos)):
    # Initiate the array of information
    array_info = np.array([saved_infos[i]])
    list_saved_things.append(array_info)
# Note how many type of variables are saved
N_saved_things = len(list_saved_things)

# Also a GUI to show that
N_pts = len(list_saved_things[0])
x_axis      = np.linspace(0, N_pts*time_btw_frame_sec, N_pts)
gui_show_saved = GUIMonitorVariables(x_axis, 'Time (sec)',
                                     list_saved_things, save_label,
                                     title='System variables')


# =============================================================================
# # Main loop
# =============================================================================
# Reference time
time_init = pygame.time.get_ticks()

# Define the following variable only to pring the time on the first iteration
# dt_elapsed = pygame.time.get_ticks() - time_init

while running:
    # A checkup
    print('(player.vx, player.vy)', (player.vx, player.vy))
    
    # Look at every event in the queue
    for event in pygame.event.get():
        # Did the user hit a key?
        if event.type == KEYDOWN:
            # Was it the Escape key? If so, stop the loop.
            if event.key == K_ESCAPE:
                running = False
        # Did the user click the window close button? If so, stop the loop.
        elif event.type == QUIT:
            running = False

    # Get the set of keys pressed and check for user input
    pressed_keys = pygame.key.get_pressed()
    
    # Update the player sprite based on user keypresses
    player.update(pressed_keys)
    
    # Update the position
    print("Warning a Coarse version of the Euler Method is Used ", pygame.time.get_ticks())
    # USe the frame duration for now
    player_x += time_btw_frame_sec*player.vx
    player_y += time_btw_frame_sec*player.vy
    #TODO Why not periodically bring back the player in a decent range, to avoid
    # infinite number for his position ?

    # Fill the screen with sky blue
    screen.fill((135, 206, 250))
    
    # Draw all sprites 
    # In this simple case only the player
    for entity in all_sprites:
        screen.blit(entity.surf, entity.rect)

    # =============================================================================
    #     # Update the portion of the map (clouds) to show in the screen
    # =============================================================================
    # In the map coordinate system, note where are the boundaries of the screen
    # Top left corner of the screen
    x_screen = player_x - SCREEN_WIDTH /2
    y_screen = player_y - SCREEN_HEIGHT/2
    # Find the index of the tile in the boundaries of the screen
    i_tile_xmin = int(round( x_screen/TILE_WIDTH ) ) - N_TILE_ADD_X
    i_tile_ymin = int(round( y_screen/TILE_HEIGHT) ) - N_TILE_ADD_Y
    # Periodic boundary
    list_ix = circularList(list_map_ix, i_tile_xmin, N_TILE_TO_SHOW_X)
    list_iy = circularList(list_map_iy, i_tile_ymin, N_TILE_TO_SHOW_Y)
    # Show the cloud in the corresponding spot
    for i in range(N_TILE_TO_SHOW_X):
        for j in range(N_TILE_TO_SHOW_Y):
            index_x = list_ix[i]
            index_y = list_iy[j]
            # BLit a cloud if there is one at this tile
            if int(map_cloud[index_x][index_y]) == int(TAG_CLOUD) :
                # Change of coord system
                pos_x = (i_tile_xmin+i)*TILE_WIDTH - x_screen
                pos_y = (i_tile_ymin+j)*TILE_HEIGHT- y_screen
                screen.blit(cloud_img, (pos_x, pos_y))
                
    # Save the data is desired
    if want_save_data:
        # Get the information
        save_label, saved_infos = player.get_info()
        for jjj in range(N_saved_things):
            prev_array = list_saved_things[jjj]
            # Initiate the array of information
            new_array = np.append(prev_array, saved_infos[jjj])
            # Check if we need to delete the history
            if len(new_array) > N_max_save:
                # Pop the first element, such that we keep only the N last
                # element (most up to date)
                new_array = np.delete(new_array, 0)
            # Overwrite the previous array with the new
            list_saved_things[jjj] = new_array
                
        # # Update the GUI to show that
        # N_pts = len(list_saved_things[0])
        # x_axis      = np.linspace(0, N_pts*time_btw_frame_sec, N_pts)
        # gui_show_saved.update(x_axis, 'Time (sec)',
        #                       list_saved_things, save_label,
        #                       title='System variables. Showing %d last data'%N_pts)

    # Update the display
    pygame.display.flip()  

    # Ensure program maintains a rate of 30 frames per second
    # Calculate the time to wait
    time_now = pygame.time.get_ticks()
    dt_elapsed = time_now - time_init
    pygame.time.wait(int(time_btw_frame - dt_elapsed))
    # Update the time for the next loop
    time_init = time_now

# Done! Time to quit.
pygame.quit()


# =============================================================================
# Plot the saved stuff
# =============================================================================

# Update the GUI to show the variables
N_pts = len(list_saved_things[0])
x_axis      = np.linspace(0, N_pts*time_btw_frame_sec, N_pts)
gui_show_saved.update(x_axis, 'Time (sec)',
                      list_saved_things, save_label,
                      title='System variables. Showing %d last data'%N_pts)

# I use this so much that it should be a function at some point
# And/or create a GUI with that

# list_y_axis = list_saved_things
# list_label  = save_label
# N_pts = len(list_saved_things[0])
# x_axis      = np.linspace(0, N_pts*time_btw_frame_sec, N_pts)


# fig, axs = plt.subplots(nrows=len(list_y_axis), 
#                         ncols=1, 
#                         sharex=True,
#                         tight_layout=True)
# for i, ax in enumerate(axs):
#     ax.plot(x_axis, list_y_axis[i], color='C%d'%i)
    
#     ax.legend([list_label[i]])
    
#     ax.set_ylabel(list_label[i])
    
# ax.set_xlabel('t (s)')
# axs[0].set_title('System variables')  

  
            
            
            
            
            
            
            
            
            
            
            
            