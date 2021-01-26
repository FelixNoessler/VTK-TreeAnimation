
import vtk

def main():
    #source
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetResolution(5)
    cylinder.SetRadius(0.2)
    cylinder.SetHeight(2)

    # mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(cylinder.GetOutputPort())

    # actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0,0.2,0.6)
    actor.RotateX(-100)
    actor.RotateY(100)

    # renderer
    render = vtk.vtkRenderer()
    render.SetBackground(1, 1, 1)
    render.AddActor(actor)

    # rendererWindow
    ren_window = vtk.vtkRenderWindow()
    ren_window.SetWindowName('Huhu')
    ren_window.SetSize(600, 600)
    ren_window.AddRenderer(render)

    # interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(ren_window)

    # initialize
    interactor.Initialize()
    ren_window.Render()
    interactor.Start()


if __name__ == '__main__':
    main()











