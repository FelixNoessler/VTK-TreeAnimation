import vtk

def main():
    # source


    # dataset 1
    points = vtk.vtkPoints()
    points.InsertNextPoint(0, 0, 0)
    points.InsertNextPoint(2, 0, 5)
    points.InsertNextPoint(2, 0, 2)
    points.InsertNextPoint(4, 0, 0)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)

    cylinder = vtk.vtkCylinderSource()
    cylinder.SetResolution(5)
    cylinder.SetRadius(0.5)
    cylinder.SetHeight(4)

    glyph3D = vtk.vtkGlyph3D()
    glyph3D.SetSourceConnection(cylinder.GetOutputPort())
    glyph3D.SetInputData(polydata)
    glyph3D.Update()

    # dataset 2
    points = vtk.vtkPoints()
    points.InsertNextPoint(0, 0, 0)
    points.InsertNextPoint(4, 0, 4)

    polydata1 = vtk.vtkPolyData()
    polydata1.SetPoints(points)

    cube = vtk.vtkCubeSource()

    glyph3D1 = vtk.vtkGlyph3D()
    glyph3D1.SetSourceConnection(cube.GetOutputPort())
    glyph3D1.SetInputData(polydata1)
    glyph3D1.Update()

    mbds = vtk.vtkMultiBlockDataSet()
    mbds.SetNumberOfBlocks(2)
    mbds.SetBlock(0, glyph3D.GetOutput())
    mbds.SetBlock(1, glyph3D1.GetOutput())



    mapper = vtk.vtkCompositePolyDataMapper2()
    mapper.SetInputDataObject(mbds)
    cdsa = vtk.vtkCompositeDataDisplayAttributes()
    mapper.SetCompositeDataDisplayAttributes(cdsa)


    mapper.SetBlockColor(1, 1,0,0)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)


    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Enable user interface interactor.
    renderer.AddActor(actor)
    renderer.SetBackground(1,1,1)
    renderWindow.Render()
    renderWindowInteractor.Start()






if __name__=='__main__':
    main()