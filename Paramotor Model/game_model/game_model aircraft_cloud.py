# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 17:37:38 2023

Goal:
    Get a strating code for the F=ma shower 

@author: mcaoue2
"""

# Just type _p() in the console 
import traceback
_p = traceback.print_last

# Import the pygame module
import pygame

# Import random for random numbers
import random

# The legendary numpy
import numpy as np

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
cloud_img= pygame.image.load("cloud.png").convert()
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
        self.surf_ori = pygame.image.load("jet.png").convert() # the .convert() call optimizes the Surface, making future .blit() calls faster.
        # The RLEACCEL constant is an optional parameter that helps pygame render more quickly on non-accelerated displays.
        # uses .set_colorkey() to indicate the color pygame will render as transparent. 
        self.surf_ori.set_colorkey((255, 255, 255), RLEACCEL)  
    
        # Center the player
        self.rect_ori = self.surf_ori.get_rect()
        self.rect_ori.centerx = (SCREEN_WIDTH -self.surf_ori.get_width() )/2
        self.rect_ori.centery = (SCREEN_HEIGHT-self.surf_ori.get_height())/2
        
        # Initiate the physical property
        self.max_v = 30 # Maximum velocity
        self.tau   = 5 # Time constant (order of magnitude of the time to reach steady state)
        # Input some initial condition
        self.v_abs = 0 # Speed, absolute
        self.v_angle = -45*np.pi/180 # Speed, angle
        self.T_angle = 0 # Direction of the trust force, also set the angle of the plane
        
        # The other physical constant are chosen to get the above property
        self.beta = 1/(self.max_v*self.tau) # A constant of the system, related to the damping
        self.mass = 0.010 # Should not affect the physic. There is a freedome in the mass, because what matters is the force/mass
        self.drag = self.mass*self.beta # It depends on the initial velocity with this drag force !
        self.max_F = self.drag*self.max_v**2 # Maximum force applied
        self.vx = self.v_abs*np.cos(self.v_angle)
        self.vy = self.v_abs*np.sin(self.v_angle)
        self.dt_solver = self.tau/100 # Increment in time for the integrator

        # Other stuff
        self.increment_angle = 5*np.pi/180
        
        # Indication of the state
        self.trust = 0 # +1, 0 or -1, for the trust force
        
        # update everything (also initiate surf and rect)
        self._update_variables()

    def model_acc(self, v, angle, fx, fy):
        # Acceleration
        cosy = np.cos(angle)
        siny = np.sin(angle)
        f_drag = -(self.drag/self.mass)*v**2*np.array([cosy, siny])
        return f_drag + np.array([fx, fy])/self.mass
        
        
    def update(self, pressed_keys):
        # The pressed keys will update some of the physical property
        # Up and down determine the direction of the force (foward or backware)
        if pressed_keys[K_DOWN]:
            self.trust = -1
        elif pressed_keys[K_UP]:
            self.trust = +1
        else:
            self.trust = 0 
            
        if pressed_keys[K_LEFT]:
            # Turn the player
            self.T_angle -= self.increment_angle
            
        if pressed_keys[K_RIGHT]:
            # Turn the player
            self.T_angle += self.increment_angle  
        
        # update everything
        self._update_variables()
        
    def _update_variables(self):
 
        # Update the physics
        # Update the acceleation based on the PREVIOUS velocity 
        # and the current force
        my_v, my_angle = get_mag_angle(self.vx, self.vy)
        fx = self.trust*self.max_F*np.cos(self.T_angle)
        fy = self.trust*self.max_F*np.sin(self.T_angle)
        acc_x, acc_y = self.model_acc(my_v, my_angle, fx, fy)
        # semi-implicit Euler method
        self.vx += acc_x * self.dt_solver
        self.vy += acc_y * self.dt_solver 
        
        # Update the image
        # Set the surface and the rect to be at a certain angle
        self.surf = pygame.transform.rotate(self.surf_ori, 
                                            -self.T_angle*rad2deg)
        self.rect = self.surf.get_rect(center = self.rect_ori.center)
                


def circularList(list_original, i_start, N_list):
    # Create a list of element from "list_original" that circulate periodically
    # from the "i_start" and N_list time. 
    # It is used to generate a periodic boundary condition on the map
    # More explanation later !
    
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
dt = 0.25 # Time interval (in whatever unit is assumed)

# Reference time
time_init = pygame.time.get_ticks()

# Main loop
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
    player_x += dt*player.vx
    player_y += dt*player.vy
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




   
            
            
            
            
            
            
            
            
            
            
            
            