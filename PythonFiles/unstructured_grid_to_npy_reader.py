import vtk
import numpy as np

from vtk.util import numpy_support

def main():

    # How many paddings should be added to file name
    PaddingAmount = 3

    Input_name, Output_name, starting_index, ending_index = get_program_parameters()

    Input_name = "../rawdata/" + Input_name
    Output_name = "../object/" + Output_name

    if (ending_index):
        for i in range(starting_index, ending_index + 1):
            file_surname = str(i).zfill(PaddingAmount)
            produce_npy(Input_name + file_surname, Output_name + file_surname)
    else:
        file_surname = str(starting_index).zfill(PaddingAmount)
        produce_npy(Input_name + file_surname, Output_name + file_surname)

def produce_npy(Input_name, Output_name):

    # Open file and read unstructured grid.
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(Input_name + ".vtk")
    reader.Update()
    grid = reader.GetOutput()
    
    # Read points.
    vtk_points = grid.GetPoints()
    print("There are ", vtk_points.GetNumberOfPoints(), " points in the data set")
    a = numpy_support.vtk_to_numpy(vtk_points.GetData())

    np.save(Output_name + ".npy", a)

def get_program_parameters():
    import argparse
    description = 'Reads an unstructured grid dataset stored in a .vtk file and output npy file which can be consumed by the rendering program later.'
    epilogue = '''
    '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('InputFile', help='E.g. experiment01/Rawdata_')
    parser.add_argument('OutputFile', help='E.g. experiment01/Object_')
    parser.add_argument('start_index', type=int, help='The starting index number, e.g. 1')
    parser.add_argument('end_index', type=int, help='The last index number, e.g. 50. If omitted only a single model will be produced', nargs='?')

    args = parser.parse_args()
    return args.InputFile, args.OutputFile, args.start_index, args.end_index

if __name__ == '__main__':
    main()