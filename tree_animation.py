import vtk
import numpy as np
import time
import math

class TreeAnimation():
    """Class that contains functions to create a vtk Animation of forest growth"""
    def __init__(self, settings, timer_set):
        """Initialize the TreeAnimation

        Creates empty lists for actors and trees.
        One tree have two actors. One for the trunk cylinder
        and another one for the crown cylinder.

        :param settings: (dictionary) settings for the forest model.
        :param timer: (dictionary) number of steps for the animation
        """

        self.set = settings
        self.timer_set = timer_set
        self.trees = []


    def initialize_render(self):
        """ Renders the tree and the grid actors """

        # vtkRenderer
        self.renderer = vtk.vtkRenderer()

        # vtkRenderWindow
        self.renderWindow = vtk.vtkRenderWindow()
        self.renderWindow.AddRenderer(self.renderer)
        self.renderWindow.SetSize(600, 600)
        self.renderWindow.SetWindowName('Forest-Modelling')

        # vtkRenderWindowInteractor
        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.renderWindow)

        # Add the grid actor
        self.renderer.AddViewProp(self.grid_actor())

        # Set the white background
        self.renderer.SetBackground(1, 1, 1)

        camera = vtk.vtkCamera();
        camera.SetPosition(self.set['gridsize_x']*1.5, 10, 10);
        #camera.SetFocalPoint(0, 0, 0);
        #self.renderer.SetActiveCamera(camera);


        # click -> red..
        style = MouseInteractorHighLightActor(self.trees)
        style.SetDefaultRenderer(self.renderer)
        self.interactor.SetInteractorStyle(style)

    def initialize_trees(self):
        self.behave = TreeBehavior(self.set, self.renderer, self.timer_set['steps'], self.trees)
        self.behave.generate_trees(self.set['number_of_trees'])


    def start_animation(self):

        # Start the timer
        timer = TimerCallback(self.timer_set['steps'], self.behave, self.renderWindow)
        self.interactor.AddObserver('TimerEvent', timer.run_timer)
        timer.timerId = self.interactor.CreateRepeatingTimer(1000)

        self.renderWindow.Render()
        self.interactor.Start()


    def grid_actor(self):
        """ Creates a 3D grid and slices the grid to a 2D grid of

        gridsize_x and gridsize_y are used to generate to grid for the forest ground.

        :return: grid_actor (vtkActor)
        """
        colors = vtk.vtkNamedColors()

        # create a grid
        xCoords = vtk.vtkFloatArray()
        for x, i in enumerate(np.linspace(0, self.set['gridsize_x']+1, self.set['gridsize_x']+1)):
            xCoords.InsertNextValue(i)

        zCoords = vtk.vtkFloatArray()
        for z, i in enumerate(np.linspace(0, self.set['gridsize_y']+1, self.set['gridsize_y']+1)):
            zCoords.InsertNextValue(i)

        # The coordinates are assigned to the rectilinear grid. Make sure that
        # the number of values in each of the XCoordinates, YCoordinates,
        # and ZCoordinates is equal to what is defined in SetDimensions().
        rgrid = vtk.vtkRectilinearGrid()
        rgrid.SetDimensions(x+1, 1, z+1)

        rgrid.SetXCoordinates(xCoords)
        rgrid.SetYCoordinates(xCoords)
        rgrid.SetZCoordinates(zCoords)

        # geometry filter to view the background grid
        plane = vtk.vtkRectilinearGridGeometryFilter()
        plane.SetInputData(rgrid)
        plane.SetExtent(0, x, 0, 0, 0, z)
        plane.Update()

        rgridMapper = vtk.vtkPolyDataMapper()
        rgridMapper.SetInputConnection(plane.GetOutputPort())

        grid_actor = vtk.vtkActor()
        grid_actor.SetMapper(rgridMapper)
        grid_actor.GetProperty().SetRepresentationToWireframe()
        grid_actor.GetProperty().SetColor(0,0,0)

        # A renderer and render window
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(colors.GetColor3d('White'))

        return grid_actor


class MouseInteractorHighLightActor(vtk.vtkInteractorStyleTrackballCamera):

    def __init__(self, trees, parent=None):
        self.trees = trees
        self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)

        self.LastPickedActor = None
        self.LastPickedProperty = vtk.vtkProperty()

    def leftButtonPressEvent(self, obj, event):

        clickPos = self.GetInteractor().GetEventPosition()

        picker = vtk.vtkPropPicker()
        picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())

        # get the new
        self.NewPickedActor = picker.GetActor()




        # If something was selected
        if self.NewPickedActor:

            x = self.NewPickedActor.GetPosition()[0]
            y = self.NewPickedActor.GetPosition()[2]


            for i,tree in enumerate(self.trees):
                x_i, y_i = tree.coord

                if x_i == x and y_i == y:
                    print(vars(self.trees[i]))
                    break


            # If we picked something before, reset its property
            if self.LastPickedActor:
                self.LastPickedActor.GetProperty().DeepCopy(self.LastPickedProperty)

            # Save the property of the picked actor so that we can
            # restore it next time
            self.LastPickedProperty.DeepCopy(self.NewPickedActor.GetProperty())
            # Highlight the picked actor by changing its properties
            self.NewPickedActor.GetProperty().SetColor([1,0,0])
            #self.NewPickedActor.GetProperty().SetDiffuse(1.0)
            #self.NewPickedActor.GetProperty().SetSpecular(0.0)
            #self.NewPickedActor.GetProperty().EdgeVisibilityOn()

            # save the last picked actor
            self.LastPickedActor = self.NewPickedActor

        self.OnLeftButtonDown()
        return


class TimerCallback():
    def __init__(self, steps, tree_behave, renWin):
        self.timer_count = 0
        self.steps = steps
        self.tree_behave = tree_behave
        self.renWin = renWin
        self.timerId = None


    def run_timer(self, obj, event):
        step = 0

        while step < self.steps:
            self.tree_behave.behave(step)

            obj.GetRenderWindow().Render()
            self.timer_count += 1
            step += 1

            if False:
                png_writer = vtk.vtkPNGWriter()
                filename = f'{step}_step.png'
                print(filename)
                png_writer.SetFileName(filename)
                a = vtk.vtkWindowToImageFilter()
                a.SetInput(self.renWin);
                png_writer.SetInputData(a.GetOutput())
                png_writer.Write()
                #stlWriter.SetInputConnection(sphereSource.GetOutputPort())
                #stlWriter.Write()


            #print(self.trees)
            #print(len(self.actor))

            #time.sleep(0.05)


        if self.timerId:
            obj.DestroyTimer(self.timerId)

class TreeBehavior:
    def __init__(self, set, renderer, steps, trees):

        self.set = set
        self.renderer = renderer
        self.steps = steps

        self.trees = trees
        self.tree_actors = []

        self.trunk_del, self.crown_del, self.tree_del = [], [], []

    def generate_trees(self, number_of_trees):
        """ Initialize trees and for each tree two vtk actors

        Trees are generated with the Tree class and filled in the tree list.
        The tree instance variables are used to generate two actors for each tree.
        The first actor is for the trunk cylinder, the second one for crown cylinder.
        """

        ### initialize trees
        before_treelenght = len(self.trees)

        for i in range(number_of_trees):
            coord = np.random.random(2)
            coord = coord[0] * self.set['gridsize_x'], \
                         coord[1] * self.set['gridsize_y']

            dist = self.calculate_dist(coord)
            treshold = 2
            great_distance = np.all(dist > treshold)
            
            if great_distance:
                tree_i = Tree(self.set, coord)
                self.trees.append(tree_i)

        ### initialize actors
        before_actorlength = len(self.tree_actors)

        for i in range(before_treelenght, len(self.trees)):
            # trunk

            r = np.random.randint(6,10)/10
            brown_color = (r, r*0.6, r*0.2)

            trunk_actor = Cylinder(self.trees[i].trunk_radius, self.trees[i].trunk_height,
                                   ( self.trees[i].coord[0],
                                     self.trees[i].trunk_height/2,
                                     self.trees[i].coord[1] ),
                                   brown_color ).cylinder_actor()
            self.tree_actors.append(trunk_actor)

            # crown


            green_color = (0, np.random.randint(3,10)/10, 0)

            crown_actor = Cylinder(self.trees[i].crown_radius, self.trees[i].crown_radius,
                                   ( self.trees[i].coord[0],
                                     self.trees[i].trunk_height,
                                     self.trees[i].coord[1] ),
                                   green_color).cylinder_actor()
            self.tree_actors.append(crown_actor)

        # Add all tree actors to the renderer
        for i in range(before_actorlength, len(self.tree_actors)):
            self.renderer.AddViewProp(self.tree_actors[i])

    def behave(self, step):
        for i in range(len(self.trees)):
            # iterator index
            trunk_i = (i + 1) * 2 - 2
            crown_i = (i + 1) * 2 - 1

            self.growth(i, trunk_i, crown_i, step)
            self.is_tree_death(trunk_i, crown_i, i)

        self.seedling()
        self.mortality()


    def growth(self, i, trunk_i, crown_i, step):
        ### growth

        dist = self.calculate_dist(self.trees[i].coord, i)
        #print(dist)

        if len(dist) > 0:
            min_dist = 2
            max_dist = min_dist * 10

            # set a growth limit, a tree can't grow better when another tree is 200m away
            next_tree = min(dist.min(), max_dist)

            # scale value from 0 to 1
            d =  (next_tree - min_dist) * 1/ (max_dist - min_dist)
        else:
            d = 1

        self.trees[i].year += 1
        scale = self.trees[i].scale + d  *4 #* self.trees[i].growth

        self.trees[i].scale = scale

        # growth of trunk
        self.tree_actors[trunk_i].SetScale(scale)
        trunk_pos = list(self.tree_actors[trunk_i].GetPosition())
        trunk_pos[1] =  self.trees[i].trunk_height * scale / 2
        self.tree_actors[trunk_i].SetPosition(trunk_pos)

        # growth of crown
        self.tree_actors[crown_i].SetScale(scale)
        crown_pos = list(self.tree_actors[crown_i].GetPosition())
        crown_pos[1] = self.trees[i].trunk_height * scale
        self.tree_actors[crown_i].SetPosition(crown_pos)




    def is_tree_death(self, trunk_i, crown_i, i):
        ### mortality
        ymax = self.tree_actors[crown_i].GetBounds()[3]

        if ymax >=  self.set['max_treesize']:
            self.tree_del.insert(0, i)
            self.trunk_del.insert(0, trunk_i)
            self.crown_del.insert(0, crown_i)

    def mortality(self):
        for index in range(0, len(self.tree_del)):
            self.renderer.RemoveViewProp(self.tree_actors[self.trunk_del[index]])
            self.renderer.RemoveViewProp(self.tree_actors[self.crown_del[index]])
            del self.trees[self.tree_del[index]]
            del self.tree_actors[self.crown_del[index]]
            del self.tree_actors[self.trunk_del[index]]

        self.trunk_del, self.crown_del, self.tree_del = [], [], []

    def seedling(self,):
        number_of_seeds = self.set['number_of_trees'] - len(self.trees)
        if number_of_seeds == 0:
            return
        else:
            self.generate_trees(number_of_seeds)


    def calculate_dist(self, coordinates, not_i = None):
        dist = np.array([])

        if not_i is not None:
            trees = self.trees[:not_i] + self.trees[not_i+1:]
        else:
            trees = self.trees

        for tree in trees:
            x = tree.coord[0]
            y = tree.coord[1]

            a = abs(x - coordinates[0])
            b = abs(y- coordinates[1])
            c =  math.sqrt(a**2 + b**2)
            dist = np.append(dist, c)

        return dist



class Cylinder():
    def __init__(self, radius, height, loc, color, opacity=1):
        self.radius = radius
        self.height = height
        self.loc = loc
        self.color = color
        self.opacity = opacity


    def cylinder_actor(self):
        """Creates a vtkCylinder, mapper and actor and set some initial settings

        The cylinder is discribed by the radius and the height.
        For each cylinder there is one mapper and one actor,
        only the actor of the cylinder is returned.
        The actor sets the location, color and opacity of the cyclinder.

        :return: actor (vtkActor)
        """

        # vtkCylinder
        cylinder = vtk.vtkCylinderSource()
        cylinder.SetResolution(100)
        cylinder.SetRadius(self.radius)
        cylinder.SetHeight(self.height)

        # mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(cylinder.GetOutputPort())

        # actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.SetPosition(self.loc)
        actor.GetProperty().SetColor(self.color)
        actor.GetProperty().SetOpacity(self.opacity)

        return actor

class Tree():
    """ Describes a tree object"""
    def __init__(self, tree_settings, coord):
        """Initialization of the tree objects

        All variables are calculated with the random number generator
        Boundaries of the variables are given in the tree_settings dictionary
        """
        self.year = 0

        # location
        self.coord = coord

        # trunk
        self.trunk_height = 1
        #self.trunk_height = max(np.random.random() * tree_settings['max_trunk_len'],
        #                      tree_settings['min_trunk_len'])

        self.trunk_radius = 0.04

        #self.trunk_radius = max(np.random.random() * tree_settings['max_trunk_rad'],
        #                     tree_settings['min_trunk_rad'])

        # crown
        self.crown_height = 0.3
        self.crown_radius = 0.2
        #self.crown_height = max(np.random.random() * tree_settings['max_crown_len'],
        #                     tree_settings['min_crown_len'])
        #self.crown_radius = max(np.random.random() * tree_settings['max_crown_rad'],
        #                     tree_settings['min_crown_rad'])

        # growth rate
        self.growth = np.random.random() * 0.7

        self.scale = 1


