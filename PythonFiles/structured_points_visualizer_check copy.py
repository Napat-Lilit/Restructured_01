import vtk
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkFiltersModeling import vtkFillHolesFilter
from vtkmodules.vtkIOLegacy import (
    vtkStructuredPointsReader
)
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkPolyDataMapper,
    vtkProperty,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)

def main():
    colors = vtkNamedColors()

    colors.SetColor('ImperialPurple', [84, 44, 93, 255])
    colors.SetColor('NavyBlue', [0, 68, 129, 255])
    colors.SetColor('GoldenRod', [224, 165, 38, 255])

    file_name, threshold = get_program_parameters()

    # Load data
    reader = vtkStructuredPointsReader()
    
    # reader = vtkUnstructuredGridReader()
    reader.SetFileName("../rawdata/" + file_name)

    mc = vtk.vtkMarchingContourFilter()
    mc.ComputeNormalsOn()

    mc.SetInputConnection(reader.GetOutputPort())
    mc.ComputeNormalsOn()
    mc.ComputeGradientsOn()
    mc.SetValue(0, threshold)  # second value acts as threshold

    # fill = vtkFillHolesFilter()
    # fill.SetInputConnection(mc.GetOutputPort())
    # fill.SetHoleSize(100.0)
    # fill.Update()

    # This has to be given, need some IO adjustment again soon
    TinyVal = 0.01
    BoundarySize = 255.0

    plane1 = vtk.vtkPlane()
    plane1.SetNormal(1, 0, 0)
    plane1.SetOrigin(TinyVal, TinyVal, TinyVal)
    plane2 = vtk.vtkPlane()
    plane2.SetNormal(0, 1, 0)
    plane2.SetOrigin(TinyVal, TinyVal, TinyVal)
    plane3 = vtk.vtkPlane()
    plane3.SetNormal(0, 0, 1)
    plane3.SetOrigin(TinyVal, TinyVal, TinyVal)
    plane4 = vtk.vtkPlane()
    plane4.SetNormal(-1, 0, 0)
    plane4.SetOrigin(BoundarySize - TinyVal, BoundarySize - TinyVal, BoundarySize - TinyVal)
    plane5 = vtk.vtkPlane()
    plane5.SetNormal(0, -1, 0)
    plane5.SetOrigin(BoundarySize - TinyVal, BoundarySize - TinyVal, BoundarySize - TinyVal)
    plane6 = vtk.vtkPlane()
    plane6.SetNormal(0, 0, -1)
    plane6.SetOrigin(BoundarySize - TinyVal, BoundarySize - TinyVal, BoundarySize - TinyVal)

    coll = vtk.vtkPlaneCollection()
    coll.AddItem(plane1)
    coll.AddItem(plane2)
    coll.AddItem(plane3)
    coll.AddItem(plane4)
    coll.AddItem(plane5)
    coll.AddItem(plane6)
    clip = vtk.vtkClipClosedSurface()
    clip.SetClippingPlanes(coll)
    clip.SetInputConnection(mc.GetOutputPort())
    clip.Update()

    # Create a mapper
    mapper = vtkPolyDataMapper()
    # mapper.SetInputConnection(mc.GetOutputPort())
    # mapper.SetInputConnection(fill.GetOutputPort())
    mapper.SetInputConnection(clip.GetOutputPort())
    mapper.ScalarVisibilityOff()

    # Visualize
    actor = vtkActor()
    actor.GetProperty().SetColor(colors.GetColor3d('ImperialPurple'))
    back_prop = vtkProperty()
    back_prop.SetDiffuseColor(colors.GetColor3d('NavyBlue'))
    actor.SetBackfaceProperty(back_prop)
    actor.SetMapper(mapper)

    renderer = vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(colors.GetColor3d('GoldenRod'))
    renderer.GetActiveCamera().SetViewUp(0.0, 0.0, 1.0)
    renderer.GetActiveCamera().SetPosition(0.0, 1.0, 0.0)
    renderer.GetActiveCamera().SetFocalPoint(0.0, 0.0, 0.0)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Azimuth(30.0)
    renderer.GetActiveCamera().Elevation(30.0)
    ren_win = vtkRenderWindow()
    ren_win.AddRenderer(renderer)
    ren_win.SetSize(640, 480)
    ren_win.SetWindowName('Model Check')

    iren = vtkRenderWindowInteractor()
    iren.SetRenderWindow(ren_win)
    ren_win.Render()
    iren.Initialize()
    iren.Start()

def get_program_parameters():
    import argparse
    description = 'Reads a structured points dataset stored in a .vtk file, constructs a 3D model, then show the resulting model on the screen'
    epilogue = '''
    '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', help='E.g. experiment01/RawResult_01.vtk.')
    parser.add_argument('threshold', type=float, help='The threshold, e.g. 0.5')
    args = parser.parse_args()
    return args.filename, args.threshold

if __name__ == '__main__':
    main()
