import vtk
from vtkmodules.vtkCommonCore import (
    VTK_VERSION_NUMBER,
    vtkVersion
)
from vtkmodules.vtkFiltersCore import (
    vtkFlyingEdges3D,
    vtkMarchingCubes
)
from vtkmodules.vtkIOGeometry import (
    vtkSTLWriter,
    vtkOBJWriter
)
from vtkmodules.vtkCommonDataModel import vtkPolyData
from vtkmodules.vtkFiltersCore import (
    vtkDecimatePro,
)
from vtkmodules.vtkIOLegacy import vtkStructuredPointsReader

def main():

    # How many paddings should be added to file name
    PaddingAmount = 3

    # vtkFlyingEdges3D was introduced in VTK >= 8.2
    use_flying_edges = vtk_version_ok(8, 2, 0)

    # Create a 3D model using flying edges or marching cubes
    if use_flying_edges:
        try:
            mc = vtkFlyingEdges3D()
        except AttributeError:
            mc = vtkMarchingCubes()
    else:
        mc = vtkMarchingCubes()

    Input_name, Output_name, threshold, target_ratio, starting_index, ending_index = get_program_parameters()

    Input_name = "../rawdata/" + Input_name
    Output_name = "../object/" + Output_name

    if (ending_index):
        for i in range(starting_index, ending_index + 1):
            file_surname = str(i).zfill(PaddingAmount)
            produce_obj(mc, Input_name + file_surname, Output_name + file_surname, threshold, target_ratio)
    else:
        file_surname = str(starting_index).zfill(PaddingAmount)
        produce_obj(mc, Input_name + file_surname, Output_name + file_surname, threshold, target_ratio)

def produce_obj(mc, Input_name, Output_name, threshold, target_ratio):
    # Load data
    reader = vtkStructuredPointsReader()
    reader.SetFileName(Input_name + ".vtk")
    mc.SetInputConnection(reader.GetOutputPort())
    mc.ComputeNormalsOn()
    mc.ComputeGradientsOn()
    mc.SetValue(0, threshold) # second value acts as threshold. first value contour number, whatever that means...

    mc.Update()
    readerMC = mc.GetOutput()

    decimate = vtkDecimatePro()
    decimate.SetInputData(readerMC)
    decimate.SetTargetReduction(target_ratio)
    decimate.PreserveTopologyOn()

    decimate.Update()
    decimated = decimate.GetOutput()

    stlWriter = vtkSTLWriter()
    stlWriter.SetFileName(Output_name + ".stl")
    stlWriter.SetInputData(decimated)
    # stlWriter.SetInputConnection(clip.GetOutputPort())
    stlWriter.Write()

def get_program_parameters():
    import argparse
    description = 'Reads a structured points dataset stored in a .vtk file and output 3D models.'
    epilogue = '''
    '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('InputFile', help='E.g. experiment01/Rawdata_')
    parser.add_argument('OutputFile', help='E.g. experiment01/Object_')
    parser.add_argument('threshold', type=float, help='The threshold, e.g. 0.8.')
    parser.add_argument('reductionratio', type=float, help='The target reduction ratio, e.g. 0.8.')
    parser.add_argument('start_index', type=int, help='The starting index number, e.g. 1')
    parser.add_argument('end_index', type=int, help='The last index number, e.g. 50. If omitted only a single model will be produced', nargs='?')

    args = parser.parse_args()
    return args.InputFile, args.OutputFile, args.threshold, args.reductionratio, args.start_index, args.end_index

# to check weather the vtk version is recent enough to use the newer flying edge or not
def vtk_version_ok(major, minor, build) :
    # return true if the requested version is at least equal to that of the actual version
    needed_version = 10000000000 * int(major) + 100000000 * int(minor) + int(build)
    try:
        vtk_version_number = VTK_VERSION_NUMBER
    except AttributeError:  # as error:
        ver = vtkVersion()
        vtk_version_number = 10000000000 * ver.GetVTKMajorVersion() + 100000000 * ver.GetVTKMinorVersion() \
                             + ver.GetVTKBuildVersion()
    if vtk_version_number >= needed_version:
        return True
    else:
        return False

if __name__ == '__main__':
    main()