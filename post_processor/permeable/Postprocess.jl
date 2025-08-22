using LinearAlgebra
using Octavian
using PyCall
using MKL
using DelimitedFiles

# Read the coordinate and element connection data from VTK file
function read_mesh(filename::String)
    # Initialize arrays
    points = []
    elements = []
    
    # Open the file for reading
    open(filename, "r") do file
        # Read through the file line by line
        while !eof(file)
            line = readline(file)
            
            # Read points
            if startswith(line, "POINTS")
                n_points = parse(Int, split(line)[2])
                while !eof(file)
                    line = strip(readline(file))
                    if isempty(line) || any(isletter, line)
                        break
                    end
                    append!(points, parse.(Float64, split(line)))
                end
                points = reshape(points, (3,:))
            end
            
            # Read elements
            if startswith(line, "CONNECTIVITY")
                while !eof(file)
                    line = strip(readline(file))
                    if isempty(line) || any(isletter, line)
                        break
                    end
                    append!(elements, parse.(Int, split(line)))
                end
                elements = reshape(elements, (4,:)) .+ 1
            end
        end
    end
    return points, elements
end

# read the v filed from VTK file
function read_v_field(filename::String)
    # Initialize displacement array
    displacement = []
    
    # Open the file for reading
    open(filename, "r") do file
        # Read through the file line by line
        while !eof(file)
            line = readline(file)
            
            # Read nodal displacement
            if startswith(line, "FIELD")
                # Skip to the line containing 'v x xxx double'
                while !startswith(line, "v ")
                    line = readline(file)
                end
                while !eof(file)
                    line = strip(readline(file))
                    if isempty(line) || (any(isletter, replace(line, "e" => "")) && !startswith(line, "v "))
                        break
                    end
                    append!(displacement, parse.(Float64, split(line)))
                end
                break
            end
        end
    end
    
    return Float64.(reshape(displacement,(1,:)))
end

# Compute the shape functions for a 2D quadrilateral element
function quad_shape_functions(xi, eta)
     
    N = [0.25 * (1 - xi) * (1 - eta),
    0.25 * (1 + xi) * (1 - eta),
    0.25 * (1 + xi) * (1 + eta),
    0.25 * (1 - xi) * (1 + eta)]

    return N'
end

# Compute the derivatives of shape functions for a 2D quadrilateral element
function quad_shape_function_derivatives(xi, eta)

    dN_dxi = [-0.25 * (1 - eta),
              0.25 * (1 - eta),
              0.25 * (1 + eta),
             -0.25 * (1 + eta)]

    dN_deta = [-0.25 * (1 - xi),
               -0.25 * (1 + xi),
                0.25 * (1 + xi),
                0.25 * (1 - xi)]

    dN = hcat(dN_dxi, dN_deta)

    return dN'
end

# Compute the Jacobian matrix for a 2D quadrilateral element at the given xi and eta.
function compute_jacobian(xi, eta, element_coords)
    dN = quad_shape_function_derivatives(xi, eta)
    dN_dxi = dN[1,:]
    dN_deta = dN[2,:]

    jacobian = zeros(Float64, 2, 2)
    for i = 1:4
        jacobian[1,1] += dN_dxi[i] * element_coords[1,i]
        jacobian[1,2] += dN_dxi[i] * element_coords[2,i]
        jacobian[2,1] += dN_deta[i] * element_coords[1,i]
        jacobian[2,2] += dN_deta[i] * element_coords[2,i]
    end

    return jacobian
end

# Compute the derivatives of shape functions in the physical space for a 2D quadrilateral element at the given xi and eta.
function shape_function_derivatives_physical(xi, eta, element_coords)

    dN = quad_shape_function_derivatives(xi, eta)
    dN_dxi = dN[1,:]
    dN_deta = dN[2,:]

    jacobian = compute_jacobian(xi, eta, element_coords)
    inv_jacobian = inv(jacobian)
    dN_dx = inv_jacobian[1,1] * dN_dxi + inv_jacobian[1,2] * dN_deta
    dN_dy = inv_jacobian[2,1] * dN_dxi + inv_jacobian[2,2] * dN_deta

    dNdx= hcat(dN_dx, dN_dy)

    return dNdx'
    
end

# Collect the shape functions
function collect_shape_functions(element, points)
    # define Gauss points
    element_coords = points[:,element]
    g = zeros(Float64, 2, 2)
    gauss_points = [-1/sqrt(3)  -1/sqrt(3);
                    1/sqrt(3)  -1/sqrt(3);
                    1/sqrt(3)   1/sqrt(3);
                   -1/sqrt(3)  1/sqrt(3)]
    gauss_weights = [1 1 1 1]
    # define the parameter coordinate to calculate the derivatives at node
    parametric_coords = [-1 -1;
                          1 -1;
                          1  1;
                         -1  1]
    # calculate the shape functions and derivatives at the Gauss points
    N = zeros(Float64,4,4)
    dN = zeros(Float64,2,4,4)
    dNdx = zeros(Float64,2,4)
    jacobian = zeros(Float64,2,2)

    for i=1:4
        N[:,i] = quad_shape_functions(gauss_points[i,1], gauss_points[i,2]) * gauss_weights[i]
        dN[:,:,i] = quad_shape_function_derivatives(gauss_points[i,1], gauss_points[i,2]) * gauss_weights[i]
        jacobian = compute_jacobian(gauss_points[i,1], gauss_points[i,2], element_coords)
        # the gradient of the shape functions in the physical space
        dNdx_all = shape_function_derivatives_physical(gauss_points[i,1], gauss_points[i,2], element_coords)
        dNdx[:,i] = dNdx_all[:,i]'
    end

    J = det(jacobian)

    return N, dN, dNdx, J
end

# compute the configural force
function configural_force!(element, N, dN, dNdx, J, v_values, dvdt_values, G_frac, G_v)
    # material parameters
    Gc = 20.2
    l = 1.0
    Lv=1/1000.0

    # Compute the fracture energy density
    v = v_values[element]
    dvdt = dvdt_values[element]

    # loop through nodes
    @inbounds @fastmath for i = 1:4
        # loop through Gauss points
        for ii = 1:4
            v_elem = N[:,ii]' * v
            dvdt_elem = N[:,ii]' * dvdt
            dv_elem = dN[:,:,ii] * v

            # calculate G_frac
            integrand = (Gc/(2*l))*(1-v_elem)^2+(Gc*l/2) * (dv_elem[1]*dv_elem[1]+dv_elem[2]*dv_elem[2])
            I_2 = Matrix{Float64}(I, 2, 2)
            psi_frac = integrand * I_2

            g_f = (psi_frac - Gc*l*(dv_elem*dv_elem'))*dN[:,i]

            G_frac[:,element[i]] =@. G_frac[:,element[i]] + g_f * J

            # calculate G_v
            g_v_gauss = @. dvdt_elem/Lv * dv_elem

            G_v[:,element[i]] = @. G_v[:,element[i]] + g_v_gauss * J

        end
    end

end

# Select nodes inside the circle and store node ID
function select_nodes_in_circle(points, elements, radius, center)
    distances = [norm(points[1:2,i] - center) for i in axes(points, 2)]
    selected_nodes_ID = findall(distances .<= radius)
    selected_nodes = points[:, selected_nodes_ID]

    # find elements that have all the nodes inside the circle
    selected_elements_ID = []
    for i in axes(elements,2)
        if all([node in selected_nodes_ID for node in elements[:,i]])
            push!(selected_elements_ID, i)
        end
    end

    return selected_nodes_ID, selected_nodes, selected_elements_ID

end

# Find the cracked nodes
function find_crack(points,v)
    # Select nodes with v>0.8
    indices = findall(v .> 0.9)
    cracked_nodes_ID = [ind[2] for ind in indices]
    cracked_nodes = points[:,cracked_nodes_ID]
    max_x_index = findmax(cracked_nodes[1,:])[2]
    node_with_max_x = cracked_nodes[:,max_x_index]
    crack_tip = node_with_max_x

    return  convert(Matrix{Float64},cracked_nodes), convert(Vector{Float64},crack_tip)
end

# Plot the mesh and selected nodes
# Initialize the Python environment
py"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.patches as patches

def plot_mesh_python(crack_tip, radius, nodal_coordinates, connectivity, selected_nodes, cracked_nodes, selected_elements_ID, time_step):
    nodal_coordinates = np.array(nodal_coordinates)
    connectivity = np.array(connectivity)-1
    selected_elements_ID = np.array(selected_elements_ID)-1
    selected_nodes = np.array(selected_nodes)
    cracked_nodes = np.array(cracked_nodes)

    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Create a list of polygons for each element
    # selected_elements = connectivity[selected_elements_ID]
    # polygons = [nodal_coordinates[element,0:2] for element in selected_elements]
    polygons = [nodal_coordinates[element,0:2] for element in connectivity]
    
    # Use PolyCollection to plot all elements at once
    poly_collection = PolyCollection(polygons, edgecolors='lightblue', linewidths=0.1, facecolors='none')
    ax.add_collection(poly_collection)
    
    # Plot the selected nodes
    ax.scatter(selected_nodes[:, 0], selected_nodes[:, 1], c='b', s=0.1, marker='o', label='Selected Nodes')

    # Plot the integration area
    circle = patches.Circle((crack_tip[0], crack_tip[1]), radius, edgecolor='blue', facecolor='none')
    ax.add_patch(circle)

    # # Plot the cracked nodes
    ax.scatter(cracked_nodes[:, 0], cracked_nodes[:, 1], c='r', s=0.1, marker='o', label='Cracked Nodes')
    
    ax.autoscale_view()
    ax.set_aspect('equal')
    ax.legend()
    plt.tight_layout()

    figname = 'mesh_with_selected_nodes' + str(time_step) + '.svg'
    plt.savefig(figname)

    return
"""
# Convert your Julia data to Python and call the function
function plot_mesh(crack_tip, radius, points, elements,selected_nodes,cracked_node, selected_elements_ID, time_step)
    py_plot_mesh = py"plot_mesh_python"
    py_plot_mesh(crack_tip, radius, points', elements', selected_nodes',cracked_node', selected_elements_ID', time_step)
end

# Usage:
function main()
    points, elements = read_mesh("T1.vtu")
    N, dN, dNdx, J = collect_shape_functions(elements[:,1], points)
    dt = 0.05

    results = []
    initial_crack = 40
    for time_step = 1:2
        # identify the cracked nodes
        v = read_v_field("T$time_step.vtu")
        v_n = read_v_field("T$(time_step-1).vtu")
        dvdt = @. (v-v_n)/dt
        cracked_node, crack_tip = find_crack(points,v)
        # select nodes for J integration
        radius = 4.9
        selected_nodes_ID, selected_nodes, selected_elements_ID = select_nodes_in_circle(points, elements, radius, crack_tip[1:2])
        # plot the mesh
        # if crack_tip[1] - initial_crack > 0.1
        #     plot_mesh(crack_tip, radius, points,elements, selected_nodes,cracked_node, selected_elements_ID, time_step)
        #     initial_crack = crack_tip[1]
        # end
        
        G_frac = zeros(Float64,2,size(points,2))
        G_v = zeros(Float64,2,size(points,2))

        for i in axes(elements,2)
            element = elements[:,i]
            configural_force!(element, N, dN, dNdx, J, v, dvdt, G_frac, G_v)
        end
        J_frac = sum(G_frac[1,selected_nodes_ID])
        J_v = sum(G_v[1,selected_nodes_ID])

        push!(results,[time_step, crack_tip[1], J_frac, J_v])
        print("time step:", time_step, "\tCrack tip:", crack_tip[1], "\tJ_frac:", J_frac, "\tJ_v:", J_v, '\n') 

        if crack_tip[1] >90
            break
        end

    end
    writedlm("configural_force.txt", results, '\t') 
end

main()