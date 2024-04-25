import vtk
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import (
    VTK_VERSION_NUMBER,
    vtkVersion
)
from vtkmodules.vtkIOLegacy import (
    vtkStructuredGridReader,
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
    reader = vtkStructuredGridReader()
    
    # reader = vtkUnstructuredGridReader()
    reader.SetFileName("../rawdata/" + file_name)

    mc = vtk.vtkMarchingContourFilter()
    mc.ComputeNormalsOn()

    mc.SetInputConnection(reader.GetOutputPort())
    mc.ComputeNormalsOn()
    mc.ComputeGradientsOn()
    mc.SetValue(0, threshold)  # second value acts as threshold

    # Create a mapper
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(mc.GetOutputPort())
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
