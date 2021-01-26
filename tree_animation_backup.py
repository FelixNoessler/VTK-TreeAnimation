import vtk
import numpy as np

class TreeAnimation():
    def __init__(self, number_of_trees, location,
                 x_len, y_len,
                 trunk_height, trunk_radius,
                 crown_height, crown_radius,
                 growth_rate,
                 timesteps):
        super().__init__()
        self.number_of_trees = number_of_trees
        self.location = location
        self.x_len = x_len
        self.y_len = y_len
        self.trunk_height = trunk_height
        self.trunk_radius = trunk_radius
        self.crown_height = crown_height
        self.crown_radius = crown_radius
        self.growth_rate = growth_rate
        self.timesteps = timesteps

    def generate_trees(self):
        self.actor = []

        for i in range(0, self.number_of_trees):

            # trunk
            trunk_actor = Cylinder(self.trunk_radius[i], self.trunk_height[i],
                                   (self.location[0,i], self.trunk_height[i]/2, self.location[1,i]),
                                   (0.5, 0.2, 0.1)).generate()
            self.actor.append(trunk_actor)

            # crown
            crown_actor = Cylinder(self.crown_radius[i], self.crown_radius[i],
                                   [self.location[0,i],self.trunk_height[i],self.location[1,i]],
                                   (0,0.9,0)).generate()
            self.actor.append(crown_actor)


    def generate_grid(self):
        colors = vtk.vtkNamedColors()

        # create a grid
        xCoords = vtk.vtkFloatArray()
        for x, i in enumerate(np.linspace(0, self.x_len, self.x_len)):
            xCoords.InsertNextValue(i)

        yCoords = vtk.vtkFloatArray()
        for y, i in enumerate(np.linspace(0, 0, 15)):
            yCoords.InsertNextValue(i)

        zCoords = vtk.vtkFloatArray()
        for z, i in enumerate(np.linspace(0, self.y_len, self.y_len)):
            zCoords.InsertNextValue(i)

        # The coordinates are assigned to the rectilinear grid. Make sure that
        # the number of values in each of the XCoordinates, YCoordinates,
        # and ZCoordinates is equal to what is defined in SetDimensions().
        rgrid = vtk.vtkRectilinearGrid()
        rgrid.SetDimensions(x + 1, y + 1, z + 1)
        rgrid.SetXCoordinates(xCoords)
        rgrid.SetYCoordinates(yCoords)
        rgrid.SetZCoordinates(zCoords)



        # geometry filter to view the background grid
        plane = vtk.vtkRectilinearGridGeometryFilter()
        plane.SetInputData(rgrid)
        plane.SetExtent(0, x, 0, 0, 0, z)
        plane.Update()

        rgridMapper = vtk.vtkPolyDataMapper()
        rgridMapper.SetInputConnection(plane.GetOutputPort())

        wireActor = vtk.vtkActor()
        wireActor.SetMapper(rgridMapper)
        wireActor.GetProperty().SetRepresentationToWireframe()
        wireActor.GetProperty().SetColor(0,0,0)


        # A renderer and render window
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(colors.GetColor3d('White'))

        return wireActor


    def render(self):
        renderer = vtk.vtkRenderer()
        renderWindow = vtk.vtkRenderWindow()
        renderWindow.AddRenderer(renderer)
        renderWindow.SetSize(600, 600)
        renderWindow.SetWindowName('Forest-Modelling')
        renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        renderWindowInteractor.SetRenderWindow(renderWindow)

        # Enable user interface interactor.
        for ac in self.actor:
            renderer.AddActor(ac)

        renderer.AddActor(self.generate_grid())
        renderer.SetBackground(1, 1, 1)


        cb = vtkTimerCallback(self.timesteps, self.actor, self.growth_rate, self.location,
                              self.trunk_height, self.trunk_radius,
                              self.crown_height, self.crown_radius,
                              renderWindowInteractor)
        renderWindowInteractor.AddObserver('TimerEvent', cb.execute)
        cb.timerId = renderWindowInteractor.CreateRepeatingTimer(self.timesteps)

        renderWindow.Render()
        renderWindowInteractor.Start()


class vtkTimerCallback():
    def __init__(self, steps, actor, growth, loc,
                 trunk_height, trunk_radius,
                 crown_height, crown_radius,
                 iren):

        self.timer_count = 0
        self.steps = steps
        self.actor = actor
        self.growth = growth
        self.loc = loc
        self.trunk_height = trunk_height
        self.trunk_radius = trunk_radius
        self.crown_height = crown_height
        self.crown_radius = crown_radius
        self.iren = iren
        self.timerId = None

    def execute(self, obj, event):
        step = 0
        while step < self.steps:
            for i in range(0, len(self.growth)):

                scale = 1 + step / self.steps * self.growth[i]

                trunk_i = (i + 1) * 2 - 2

                self.actor[trunk_i].SetScale(scale, scale, scale)
                trunk_ymin = self.actor[trunk_i].GetBounds()[2]

                crown_i = (i + 1) * 2 - 1
                self.actor[crown_i].SetScale(scale, scale, scale)


                if trunk_ymin < 0:
                    trunk_pos = self.actor[trunk_i].GetPosition()
                    trunk_pos -= np.array([0, trunk_ymin, 0])
                    self.actor[trunk_i].SetPosition(trunk_pos)

                    crown_pos = self.actor[crown_i].GetPosition()
                    crown_pos -= np.array([0, 2*trunk_ymin, 0])
                    self.actor[crown_i].SetPosition(crown_pos)






            iren = obj
            iren.GetRenderWindow().Render()
            self.timer_count += 1
            step += 1
        if self.timerId:
            iren.DestroyTimer(self.timerId)


class Cylinder():
    def __init__(self, radius, height, loc, color, opacity=1):
        self.radius = radius
        self.height = height
        self.loc = loc
        self.color = color
        self.opacity = opacity


    def generate(self):
        cylinder = vtk.vtkCylinderSource()
        cylinder.SetResolution(50)
        cylinder.SetRadius(self.radius)
        cylinder.SetHeight(self.height)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(cylinder.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.SetPosition(self.loc)
        actor.GetProperty().SetColor(self.color)
        actor.GetProperty().SetOpacity(self.opacity)
        return actor

class Tree():
    def __init__(self, gridsize_x, grid_size_y,
                 min_trunk_len, max_trunk_len,
                 min_trunk_rad, max_trunk_rad,
                 min_crown_len, max_crown_len,
                 min_crown_rad, max_crown_rad):

        self.coord = np.random.random(2)
        self.coord = [self.coord[0] * gridsize_x, self.coord[1] * gridsize_y ]

        self.trunk_len = max(np.random.random() * max_trunk_len,
                              min_trunk_len)

        self.trunk_rad = max(np.random.random() * max_trunk_rad,
                             min_trunk_rad)

        self.crown_len = max(np.random.random() * max_crown_len,
                             min_crown_len)

        self.crown_rad = max(np.random.random() * max_crown_rad,
                             min_crown_rad)

        self.growth = np.random.random() *0.7
