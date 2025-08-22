using PyCall

py"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.patches as patches
import vtk
from vtkmodules.util import numpy_support as ns

# generate uniform grid
def create_grid(num_nodes_x=200, num_nodes_y=200):
    # Define nodes
    x = np.linspace(0, 100, num_nodes_x)
    y = np.linspace(0, 100, num_nodes_y)
    X, Y = np.meshgrid(x, y)
    nodes = np.column_stack((X.ravel(), Y.ravel()))

    # Define elements (quadrilaterals)
    elements = []
    for i in range(num_nodes_x - 1):
        for j in range(num_nodes_y - 1):
            n0 = i * num_nodes_y + j
            n1 = i * num_nodes_y + (j+1)
            n2 = (i+1) * num_nodes_y + (j+1)
            n3 = (i+1) * num_nodes_y + j
            elements.append([n0, n1, n2, n3])

    # Create VTK mesh
    points = vtk.vtkPoints()
    for node in nodes:
        points.InsertNextPoint(node[0], node[1], 0.0)

    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)

    # Add cells to the grid
    for elem in elements:
        id_list = vtk.vtkIdList()
        for nid in elem:
            id_list.InsertNextId(nid)
        grid.InsertNextCell(vtk.VTK_QUAD, id_list)

    return grid, nodes, elements

def interpolate_vtk(source_file, target_grid):
    # Load source VTU data
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(source_file)
    reader.Update()

    # Resample source data onto the generated target mesh
    resample = vtk.vtkResampleWithDataSet()
    resample.SetSourceData(reader.GetOutput())
    resample.SetInputData(target_grid)
    resample.SetPassPointArrays(True)
    resample.Update()

    # Extract the interpolated data
    resampled_vtk = resample.GetOutput()
    resampled_data_vtk = resampled_vtk.GetPointData().GetArray("v")

    # Convert VTK array to Python list
    resampled_values = [resampled_data_vtk.GetValue(i) for i in range(resampled_data_vtk.GetNumberOfTuples())]

    return resampled_values

def collect_point_data_on_line_python(vtu_file_path):
    # Step 1: Read the VTU file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_file_path)
    reader.Update()

    # Step 2: Create a line source from (0,100) to (100, 100)
    lineSource = vtk.vtkLineSource()
    lineSource.SetPoint1(0, 55, 0)
    lineSource.SetPoint2(100, 55, 0)
    lineSource.SetResolution(99)
    lineSource.Update()

    # Step 3: Probe the dataset
    probeFilter = vtk.vtkProbeFilter()
    probeFilter.SetSourceData(reader.GetOutput())
    probeFilter.SetInputConnection(lineSource.GetOutputPort())
    probeFilter.Update()

    # Step 4: Extract the point data
    probedData = probeFilter.GetOutput()
    strainArray = probedData.GetPointData().GetArray("strain")
    PArray = probedData.GetPointData().GetArray("P")

    # Convert to a NumPy array
    strain = ns.vtk_to_numpy(strainArray)
    P = ns.vtk_to_numpy(PArray)

    return strain, P

def plot_mesh_python(crack_tip, radius, nodal_coordinates, connectivity, selected_nodes, cracked_nodes, selected_elements_ID, time_step):
    nodal_coordinates = np.array(nodal_coordinates)
    connectivity = np.array(connectivity)-1
    selected_elements_ID = np.array(selected_elements_ID)-1
    selected_nodes = np.array(selected_nodes)
    cracked_nodes = np.array(cracked_nodes)

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.axis('off')
    
    # Create a list of polygons for each element
    selected_elements = connectivity[selected_elements_ID]
    polygons1 = [nodal_coordinates[element,0:2] for element in selected_elements[0]]
    polygons2 = [nodal_coordinates[element,0:2] for element in connectivity]
    
    # Use PolyCollection to plot all elements at once
    poly_collection1 = PolyCollection(polygons1, edgecolors='r', linewidths=1, facecolors='none', linestyle='dashed')
    poly_collection2 = PolyCollection(polygons2, edgecolors='b', linewidths=0.1, facecolors='none')
    ax.add_collection(poly_collection2)
    ax.add_collection(poly_collection1)
    
    # Plot the selected nodes
    # ax.scatter(selected_nodes[:, 0], selected_nodes[:, 1], c='r', s=5, marker='o', label='Selected Nodes')

    # Plot the integration area
    circle = patches.Circle((crack_tip[0], crack_tip[1]), radius, edgecolor='r', facecolor='none')
    ax.add_patch(circle)

    # # Plot the cracked nodes
    # ax.scatter(cracked_nodes[:, 0], cracked_nodes[:, 1], c='b', s=2, marker='o', label='Cracked Nodes')
    
    ax.autoscale_view()
    ax.set_aspect('equal')
    # ax.legend()
    plt.tight_layout()

    figname = 'mesh_with_selected_nodes' + str(time_step) + '.svg'
    plt.savefig(figname)

    plt.close(fig)

    return

"""