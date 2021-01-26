import vtk

counter = 0
max_count = 6


def main():
    sphere = vtk.vtkSphereSource()
    sphere.Update()

    filter = vtk.vtkProgrammableFilter()
    filter.SetInputConnection(sphere.GetOutputPort())

    filter.SetExecuteMethod(adjust_points(filter))

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(filter.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    renderer = vtk.vtkRenderer()

    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetWindowName("DataAnimation")

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)

    interactor.Initialize()
    interactor.CreateRepeatingTimer(500)

    #interactor.AddObserver('TimerEvent', filter)

    renderer.AddActor(actor)
    render_window.Render()
    interactor.Start()



def adjust_points(filter):

    inPts = filter.GetPolyDataInput().GetPoints();
    numPts = inPts.GetNumberOfPoints();

    new_points = vtk.vtkPoints()
    new_points.SetNumberOfPoints(numPts)

    for i in range(numPts):
        pass

if __name__ == '__main__':
    main()