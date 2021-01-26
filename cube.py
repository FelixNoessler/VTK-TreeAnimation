import vtk

# source
cube = vtk.vtkCubeSource()
cube.SetXLength(2)
cube.SetYLength(1)
cube.SetZLength(0.5)

# no filter
# mapper
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(cube.GetOutputPort())

#actor
actor = vtk.vtkActor()
actor.SetMapper(mapper)

# renderer
render = vtk.vtkRenderer()
render.SetBackground(0,0,0)
render.AddActor(actor)

# rendererWindow
ren_window = vtk.vtkRenderWindow()
ren_window.SetWindowName('Huhu')
ren_window.SetSize(600,600)
ren_window.AddRenderer(render)

# interactor
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(ren_window)

# initialize
interactor.Initialize()
ren_window.Render()
interactor.Start()