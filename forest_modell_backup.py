import numpy as np
import tree_animation

number_of_trees = 100

# location of tree
gridsize_x, gridsize_y = 40,40

coordinates  =np.random.random( (number_of_trees, 2))
coordinates = [coordinates[:,0] * gridsize_x, coordinates[:,1] * gridsize_y ]
coordinates = np.array(coordinates)



# size of tree
min_trunk_len, max_trunk_len = 3,10
trunk_len = np.random.random(number_of_trees) * max_trunk_len
trunk_len[trunk_len < min_trunk_len] = min_trunk_len
trunk_len = np.sort(trunk_len)

min_trunk_rad, max_trunk_rad = 0.2, 0.4
trunk_rad = np.random.random(number_of_trees) * max_trunk_rad
trunk_rad[trunk_rad < min_trunk_rad] = min_trunk_rad
trunk_rad = np.sort(trunk_rad)

min_crown_len, max_crown_len = 1,5
crown_len = np.random.random(number_of_trees)  * max_crown_len
crown_len[crown_len < min_crown_len]  = min_crown_len
crown_len = np.sort(crown_len)

min_crown_rad, max_crown_rad = 0.5, 3
crown_rad = np.random.random(number_of_trees) * max_crown_rad
crown_rad[crown_rad < min_crown_rad] = min_crown_rad
crown_rad = np.sort(crown_rad)


growth_rate = np.random.random(number_of_trees) * 0.7

t = tree_animation.TreeAnimation(number_of_trees, coordinates,
                                 gridsize_x, gridsize_y,
                                 trunk_len, trunk_rad,
                                crown_len, crown_rad,
                                 growth_rate,
                                 400)
t.generate_trees()
t.initialize_render()