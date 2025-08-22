# Calculate from the original mesh
using LightXML
using LinearAlgebra
using Statistics
using Octavian
using DelimitedFiles

include("Pymodule.jl")
include("Input_Ferro.jl")

# get data from vtu file
function read_vtu(file_path)
    # Load the XML document
    xdoc = parse_file(file_path)
    xroot = root(xdoc)

    # Navigate to the Piece element
    unstructured_grid = find_element(xroot, "UnstructuredGrid")
    piece = find_element(unstructured_grid, "Piece")

    # Extract points
    function split_content(element)
        # Helper function to split and convert string data to float
        return parse.(Float64, split(content(element)))
    end
    points_element = find_element(piece, "Points")
    data_array = find_element(points_element, "DataArray")
    points_data_raw = split_content(data_array)
    points_data = reshape(points_data_raw, (3, div(length(points_data_raw), 3)))

    # Extract cells
    # Custom function to get child elements by tag
    function get_elements_by_tag(parent, tag)
        return [child for child in child_elements(parent) if name(child) == tag]
    end
    function split_content_int(element)
        # Helper function to split and convert string data to Int64
        return parse.(Int64, split(content(element)))
    end
    cells_element = find_element(piece, "Cells")
    data_arrays = get_elements_by_tag(cells_element, "DataArray")
    connectivity = first(filter(x -> attribute(x, "Name") == "connectivity", data_arrays))
    connectivity_data_raw = split_content_int(connectivity)
    connectivity_data = reshape(connectivity_data_raw, (4, div(length(connectivity_data_raw), 4))) .+ 1  

    # Extract point data
    point_data_element = find_element(piece, "PointData")
    scalar_data_arrays = get_elements_by_tag(point_data_element, "DataArray")
    point_data = Dict()
    for data_array in scalar_data_arrays
        name = attribute(data_array, "Name")
        point_data[name] = split_content(data_array)
    end

    return points_data, connectivity_data, point_data
end

# generate the uniform grid
function create_grid(num_nodes_x=100, num_nodes_y=100)
    py_create_grid = py"create_grid"
    grid, nodes, elements = py_create_grid(num_nodes_x, num_nodes_y)
    elements = elements .+ 1
    return grid, nodes', elements'
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
    dNdx = zeros(Float64,2,4,4)
    jacobian = zeros(Float64,2,2)

    for i=1:4
        N[:,i] = quad_shape_functions(gauss_points[i,1], gauss_points[i,2]) * gauss_weights[i]
        dN[:,:,i] = quad_shape_function_derivatives(gauss_points[i,1], gauss_points[i,2]) * gauss_weights[i]
        jacobian = compute_jacobian(gauss_points[i,1], gauss_points[i,2], element_coords)
        # the gradient of the shape functions in the physical space
        dNdx_all = shape_function_derivatives_physical(gauss_points[i,1], gauss_points[i,2], element_coords)
        # print(dNdx_all)
        dNdx[:,:,i] = dNdx_all
    end

    J = det(jacobian)

    return N, dN, dNdx, J
end

# compute the configural force
function configural_force!(element, N, dN, dNdx, J, v_values, dvdt_values, G_frac, G_v)
    # material parameters
    _, _, Fract, _ = Input_Ferro()
    Gc = Fract.cenerg
    l = Fract.constl
    Lv = Fract.Lv

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
            integrand = (Gc/(2.0*l))*v_elem^2+(Gc*l/2.0) * (dv_elem[1]*dv_elem[1]+dv_elem[2]*dv_elem[2])
            I_2 = Matrix{Float64}(I, 2, 2)
            psi_frac = integrand * I_2

            g_f = 6*(psi_frac - Gc*l*(dv_elem*dv_elem'))*dN[:,i,ii]

            G_frac[:,element[i]] = @. G_frac[:,element[i]] + g_f * J

            # calculate G_v
            g_v_gauss = @. dvdt_elem/Lv * dv_elem

            G_v[:,element[i]] = @. G_v[:,element[i]] + g_v_gauss * N[i,ii] * J

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
    cracked_nodes_ID = findall(v .> 0.9)
    cracked_nodes = points[:,cracked_nodes_ID]
    max_x_index = findmax(cracked_nodes[1,:])[2]
    node_with_max_x = cracked_nodes[:,max_x_index]
    crack_tip = node_with_max_x

    return  convert(Matrix{Float64},cracked_nodes), convert(Vector{Float64},crack_tip)
end

# resample v filed
function resample_v(source_file, grid)
    py_resample_v = py"interpolate_vtk"
    v = py_resample_v(source_file, grid)
    return v
end

# Plot the final mesh
function plot_mesh(crack_tip, radius, points, elements,selected_nodes,cracked_node, selected_elements_ID, time_step)
    py_plot_mesh = py"plot_mesh_python"
    py_plot_mesh(crack_tip, radius, points', elements', selected_nodes',cracked_node', selected_elements_ID', time_step)
end

# Calculate the load force
function load_force(vtu_file_path)
    py_collect_point_data_on_line = py"collect_point_data_on_line_python"
    strain, P = py_collect_point_data_on_line(vtu_file_path)

    stress_Y = (strain[:,2] .- 0.05 .* P[:,2] .^2 .+ 0.015 .* P[:,1] .^2) .* 1766 + (strain[:,1] .- 0.05 .* P[:,1] .^2 .+ 0.015 .* P[:,2] .^2) .* 802
    F_Y=mean(stress_Y)*100
    return F_Y
end

function G_frac!(results, time_step, initial_crack, grid, points_r, elements_r, N,dN,dNdx,J)
    points, elements, point_data = read_vtu("T$time_step.vtu")
        v = point_data["v"]
        v_r = resample_v("T$(time_step).vtu", grid)
        vn_r = resample_v("T$(time_step-1).vtu", grid)

        # identify the cracked nodes
        dt = 0.05
        dvdt = @. (v_r-vn_r)/dt
        cracked_node, crack_tip = find_crack(points,v)
        # select nodes for J integration
        radius = 5
        selected_nodes_ID, selected_nodes, selected_elements_ID = select_nodes_in_circle(points_r, elements_r, radius, crack_tip[1:2])
        # plot the mesh
        # if crack_tip[1] - initial_crack > 0.1 || time_step == 1
        #     # plot_mesh(crack_tip, radius, points_r,elements_r, selected_nodes,cracked_node, selected_elements_ID, time_step)
        #     plot_J(time_step, initial_crack)
        # end
        
        G_frac = zeros(Float64,2,size(points_r,2))
        G_v = zeros(Float64,2,size(points_r,2))

        for i in axes(elements_r,2)
            element_r = elements_r[:,i]
            configural_force!(element_r, N, dN, dNdx, J, v_r, dvdt, G_frac, G_v)
        end
        J_frac = sum(G_frac[1,selected_nodes_ID])
        J_v = sum(G_v[1,selected_nodes_ID])

        results[time_step,1:4]=[time_step, crack_tip[1], J_frac, J_v]
        
        return crack_tip

end

function plot_J(time_step, initial_crack)
    points, elements, point_data = read_vtu("T$time_step.vtu")
    v = point_data["v"]
    dt = 0.05

    cracked_node, crack_tip = find_crack(points,v)
    # select nodes for J integration
    radius = 5
    selected_nodes_ID, selected_nodes, selected_elements_ID = select_nodes_in_circle(points, elements, radius, crack_tip[1:2])
    # plot the mesh
    if crack_tip[1] - initial_crack > 0.1 || time_step == 1
        plot_mesh(crack_tip, radius, points,elements, selected_nodes,cracked_node, selected_elements_ID, time_step)
        initial_crack = crack_tip[1]
    end
    
    return crack_tip

end

# Usage:
function main()

    results = zeros(Float64, 4000, 5)
    initial_crack = 40
    grid, points_r, elements_r=create_grid(201, 201)
    N, dN, dNdx, J = collect_shape_functions(elements_r[:,1], points_r)

    useful_step = 0
    for time_step = 1:4000

        crack_tip = G_frac!(results, time_step, initial_crack, grid, points_r, elements_r, N,dN,dNdx,J)
        results[time_step,5] = load_force("T$time_step.vtu")
        print("time step:", time_step, "\tCrack tip:", crack_tip[1], "\tJ_frac:", results[time_step,3], "\tJ_v:", results[time_step,4], "\tF:", results[time_step,5], '\n') 
        # print("time step:", time_step, '\n') 

        initial_crack = crack_tip[1]
        # if crack_tip[1] >90
        #     break
        # end

        useful_step = time_step

    end
    writedlm("configural_force_original.txt", results[1:useful_step,:], '\t') 
end


main()
