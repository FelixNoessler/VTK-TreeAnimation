import numpy as np
import tree_animation

settings = {
    'number_of_trees':  50,
    'gridsize_x':       20,
    'gridsize_y':       20,
    'min_trunk_len':    2,
    'max_trunk_len':    4,
    'min_trunk_rad':    0.1,
    'max_trunk_rad':    0.2,
    'min_crown_len':    1,
    'max_crown_len':    1.3,
    'min_crown_rad':    0.5,
    'max_crown_rad':    0.7,
    'max_treesize':     10
}

timer = { 'steps' : 50 }

t = tree_animation.TreeAnimation(settings, timer)

t.initialize_render()
t.initialize_trees()
t.start_animation()