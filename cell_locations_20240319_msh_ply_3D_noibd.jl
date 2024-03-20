# using FileIO
# using Meshes
#using GLMakie
# using Printf
using LinearAlgebra
#using ColorSchemes
#import ColorSchemes.rainbow
# using WriteVTK
# using BenchmarkTools
include("build_edges_from_msh_3D.jl")
struct BC
    ID::Int64
    type::Int64
    value::Float64
end

function uniquenode(filename,cell_number,BCarray,ICvalue,ℾ,Sourceᵤ,Sourceₚ)


    npc = 3

    ## Read ply file for given cell
    num_nodes = 0
    num_cells = 0
    num_zones = 0
    num_edges = 0
    num_dim = 3
    count = 0
    # Parse header
    for i in eachline(filename)
        if startswith(i,"element vertex ")
            num_nodes = parse(Int64,chop(i, head = 15,tail = 0))
        end
        if startswith(i,"element face ")
             num_cells = parse(Int64,chop(i, head = 13,tail = 0))
        end
         count += 1
        if startswith(i,"end_header")
            break
        end
    end
    
    @show num_nodes
    @show num_cells
  

    
    # Read nodes (vertex) and cells (faces)
    nodes = Array{Float64,2}(undef, num_nodes,3)
   
    cells = Array{Int64,2}(undef, num_cells,npc*3+1)
    cellfaces = Array{Int64,1}(undef,3)
    count2 = 0
    for i in eachline(filename)
         count2 += 1
        # Nodes are stored first in the vertex section of file
        if count2 > count && count2 < count + num_nodes +1
            nodes[count2-count,1:3] = parse.(Float64,split(i))[1:3]
        end
        # Cells are stored next in the face section of the file
        if count2 > count + num_nodes
            celltemp = parse.(Int64,split(i))
            cellname = celltemp[1]
            @views cellfaces = celltemp[2:4].+1
            cellfaces = sort(cellfaces)
            @views cells[count2-count - num_nodes,1:3] .=  [cellfaces[1],cellfaces[2],cellfaces[3]]
        end
    end
 
    if npc == 3
        circle_path = [1 2;2 3;3 1]
    else
        if npc == 4
            circle_path = [1 2;2 3;3 4;4 1]
        else
            println("npc (number of cells per node) must be 3 or 4)")
            return
        end
    end

    cell = Array{Int64,2}(undef,1,npc*3+1)
    edge_to_cell = Array{Int64,2}(undef,num_cells*npc,5)
    edge_to_cell_sort = Array{Int64,2}(undef,num_cells*npc,5)
   
    for (ic,cell) in enumerate(eachrow(cells))
        
        for j in range(1,npc)
            
            if cell[circle_path[j,1]] < cell[circle_path[j,2]]
                @views  edge_to_cell[(ic-1)*npc+j,1:5] .=  [cell[circle_path[j,1]],cell[circle_path[j,2]],ic,j,0]
            else
                @views  edge_to_cell[(ic-1)*npc+j,1:5] .=  [cell[circle_path[j,2]],cell[circle_path[j,1]],ic,j,1]
            end
            
        end

        
    end
    # Sort by nodes so that rows of edge_to_cell that share a edge are adjacent
    # edge_to_cell_sort = similar(edge_to_cell)

    
    edge_to_cell_sort = sortslices(edge_to_cell,dims=1,by=x->(x[1],x[2]))
    

    node_list = Array{Int64,2}(undef,2,5)
    count1 = 0
    count2 = 0
    skip = 1
    edges = Array{Int64,2}(undef,num_cells*npc,11)
    boundary_edges = Array{Int64,2}(undef,num_cells*npc,11)
    num_edges,num_boundary_edges = edge_comparison(edge_to_cell_sort,edges,boundary_edges,cells,node_list,npc,count1,count2,skip)
    edges = edges[1:num_edges,:]
    boundary_edges = boundary_edges[1:num_boundary_edges,:]
    
  

    

    
    vertices = Array{Float64,3}(undef,num_cells,npc,3)
    t = Array{Float64,3}(undef,num_cells,npc,3)
    t̂ = Array{Float64,3}(undef,num_cells,npc,3)

    face_mid = Array{Float64,3}(undef,num_cells,npc,3)
    n_cell = Array{Float64,2}(undef,num_cells,3)
    n̂_cell = Array{Float64,2}(undef,num_cells,3)
    
    n̂_face = Array{Float64,3}(undef,num_cells,npc,3)
    
    volume_cell = Vector{Float64}(undef,num_cells)
    d_cell_face = Array{Float64,2}(undef,num_cells,npc)
    cell_center = Array{Float64,2}(undef,num_cells,3)
    

    # Each cell has some associated data.  These are:
    # The center of the cell (cell_center)
    # The vertices of the cell (vertices)
    # The tangent to each face (t), which are not unit vectors
    # The midpoints of the faces of the cell (face_mid)
    # The unit normal vector to each face (n̂_face), which are unit vectors
    # The normal vector of the cell (n_cell), a unit vector
    # The unitnormal vector of the cell (n̂_cell), a unit vector
    # The distance from cell center to each face midpoint (d_cell_face)
    # Cell area, also called cell volume in 2D, even though it is an area, is cell_volume
    # In contrast, the length of each 'face' is called an area.
    for (ic,cell) in enumerate(eachrow(cells))
        for j in range(1,npc)
            @views  vertices[ic,j,1:3] .= nodes[cell[j],1:3]
        end
    #     @show cell
        vertex_temp = vertices[ic,:,1:3]
         cell_center[ic,1:3] = sum(vertex_temp,dims=1)/npc


        for j in range(1,npc)
            @views  t[ic,j,1:3] .= vertices[ic,circle_path[j,2],:]-vertices[ic,circle_path[j,1],:]
            @views  face_mid[ic,j,1:3] .= (vertices[ic,circle_path[j,2],:]+vertices[ic,circle_path[j,1],:])/2
        end
       
       
     
        
        @views  n_cell[ic,1:3] .= cross(t[ic,1,1:3],-t[ic,2,1:3])
        @views  n̂_cell[ic,1:3] .= n_cell[ic,1:3]./norm(n_cell[ic,1:3])
        
        for j in range(1,npc)
         
           
            @views  t̂[ic,j,1:3] .= t[ic,j,1:3]./norm(t[ic,j,1:3])
            @views  n̂_face[ic,j,1:3] .= cross(t̂[ic,j,1:3],n̂_cell[ic,1:3])
           
            @views  d_cell_face[ic,j] = norm([cell_center[ic,1:3]-face_mid[ic,j,1:3]])
            
        end
       
        @views volume_cell[ic] = norm(n_cell[ic,1:3])/2 # Cell area, aka cell volume in 2D
        
        
        
    end
   
    num_bfaces = sum(edges[:,6])
    @show num_edges,num_boundary_edges,num_bfaces
    
    area_face = Vector{Float64}(undef,num_edges)
    
   
    area_bface = Vector{Float64}(undef,num_bfaces)
    lf = Array{Float64,3}(undef,num_edges,3,2)
    δ_face = Vector{Float64}(undef,num_edges)
    weight_face = Array{Float64,2}(undef,num_edges,2)
    lbf = Array{Float64,2}(undef,num_bfaces,3)
    δ_bface = Vector{Float64}(undef,num_bfaces)
    weight_bface = Vector{Float64}(undef,num_bfaces)
    for (ie,edge) in enumerate(eachrow(edges))
        node1 = edge[1]
        node2 = edge[2]
        cell1 = edge[3]
        cell2 = edge[4]
        face_cell1 = edge[7]
        face_cell2 = edge[8]
        bface_number = edge[9]
        @views area_face[ie] = norm(nodes[node2,1:3]-nodes[node1,1:3])
        if edge[6] == 0
            # weight_face depends on the cell order in edges.  In edges, the first cell is the one
            # with lower number.
            weight_face[ie,1] = d_cell_face[cell2,face_cell2]/(d_cell_face[cell1,face_cell1] + d_cell_face[cell2,face_cell2])
            weight_face[ie,2] = 1.0 - weight_face[ie,1]
            # lf also depends on the cell order in edges
            @views  lf[ie,1:3,1] .=  cell_center[cell2,1:3] .- cell_center[cell1,1:3]
            @views  lf[ie,1:3,2] .= cell_center[cell1,1:3] .- cell_center[cell2,1:3]
            # δ_face is the same for both cells
            @views  δ_face[ie] = abs(dot(lf[ie,1:3,1],n̂_face[cell1,face_cell1,1:3]))
        else
            # Only one cell per edge at boundary
            area_bface[bface_number] = area_face[ie]
            weight_face[ie,1] = 1.0
            weight_bface[bface_number,1] = 1.0
            # @show cell1,face_cell1,edge
            @views  lf[ie,1:3,1] .= face_mid[cell1,face_cell1,1:3] .- cell_center[cell1,1:3]
            @views  lbf[bface_number,1:3] .= lf[ie,1:3,1]
            @views  δ_face[ie] = abs(dot(lf[ie,1:3,1],n̂_face[cell1,face_cell1,1:3]))
            δ_bface[bface_number] = δ_face[ie]

    

        end
    end
    
   

    weight_node = Array{Float64,2}(undef,num_cells,npc)
    weight_node_sum = zeros(num_nodes)
    for (ic,cell) in enumerate(eachrow(cells))
        for j in range(1,npc)
            @views  weight_node[ic,j] = 1/norm(nodes[cell[j],1:3]-cell_center[ic,1:3])
            weight_node_sum[cell[j]] += weight_node[ic,j]
        end
    end

    for (ic,cell) in enumerate(eachrow(cells))
        for j in range(1,npc)
            weight_node[ic,j] = weight_node[ic,j]/weight_node_sum[cell[j]]
        end
    end
    

    
    phi = Vector{Float64}(undef, num_cells)


    IC(cells,phi,ICvalue)

    phi_bface = zeros(num_bfaces)
    phi_node = zeros(num_nodes)
    t̂_dot_lf = Array{Float64,2}(undef,num_cells,npc)
    aₚ = Vector{Float64}(undef,num_cells)
    aₙ = Array{Float64}(undef,num_cells,npc)
    sourceᵤ = Vector{Float64}(undef,num_cells)
    sourceₚ = Vector{Float64}(undef,num_cells)
    sourceₛ = Vector{Float64}(undef,num_cells)
    BC(boundary_edges,phi_bface,phi_node,BCarray)
    
   for (ic,cell) in enumerate(eachrow(cells))
        aₚ[ic] = 0.0
        sourceᵤ[ic] = Sourceᵤ*volume_cell[ic]
        sourceₚ[ic] = Sourceₚ*volume_cell[ic]
        for j in range(1,npc)   
            edge = cell[npc+j] 
            # @show ic,edge,cell,npc,j
            if edges[edge,6] == 0  
                aₙ[ic,j] = ℾ*area_face[edge]/δ_face[edge]
                aₚ[ic] += aₙ[ic,j]
                @views t̂_dot_lf[ic,j] = dot(t̂[ic,j,1:3],lf[edge,1:3,cell[npc*3+1]])
            else
                for bvalue in BCarray
                    aₙ[ic,j] = 0.0
                    @views t̂_dot_lf[ic,j] = dot(t̂[ic,j,1:3],lf[edge,1:3,1])
                    if edges[edge,5] == bvalue.ID
                        if bvalue.type == 2 # Neumann
                            sourceᵤ[ic] += bvalue.value*area_face[edge]
                        else
                            sourceᵤ[ic] += bvalue.value*ℾ*area_face[edge]/δ_face[edge]
                            sourceₚ[ic] += -ℾ*area_face[edge]/δ_face[edge]
                        end
                    end
                end
            end
            
        end
        aₚ[ic] -= sourceₚ[ic]
    end


    for iter in range(1,1000)
       # @show iter
        for (ic,cell) in enumerate(eachrow(cells))
            
            for j in range(1,npc)
                 phi_node[cell[j]] += phi[ic]*weight_node[ic,j]               
            end
        end
    

        # Fix Dirichlet boundary nodes
        for (i,bface) in enumerate(eachrow(boundary_edges))
            #@show i,bface
            for bvalue in BCarray
            # @show bface[5], bvalue.ID
                if bface[5] == bvalue.ID
                    # @show bvalue.type
                    if bvalue.type == 1 # Dirichlet
                        # @show bvalue.value
                         phi_node[bface[1]] = bvalue.value
                         phi_node[bface[2]] = bvalue.value
                    end
                end
            end
        end

        # Compute skew sources
        for (ic,cell) in enumerate(eachrow(cells))
            sourceₛ[ic] = 0.0
            for j in range(1,npc)
                 sourceₛ[ic] += t̂_dot_lf[ic,j]*ℾ/δ_face[cell[j+npc]]*(phi_node[cell[circle_path[j,2]]] - phi_node[cell[circle_path[j,1]]])
                # @show aₙ[ic,j],area_face[cell[npc+j]],δ_face[cell[npc+j]]
            end
            
        #    @show sourceₛ[ic],aₚ[ic],sourceᵤ[ic],sourceₚ[ic],t̂_dot_lf[ic,1:npc],ℾ
        end

        # Update solution using GS formula
        for (ic,cell) in enumerate(eachrow(cells))
            
            phi[ic] = (sourceᵤ[ic] + sourceₛ[ic])/aₚ[ic]
            for j in range(1,npc)
                
                if cell[j+npc*2] > 0
                     phi[ic] += aₙ[ic,j]*phi[cell[j+npc*2]]/aₚ[ic]
                    #@show iter,cell[j+npc*2],phi[cell[j+npc*2]],aₙ[ic,j],phi[ic],aₙ[ic,j]*phi[cell[j+npc*2]]
                end
                
                 phi_node[cell[j]] = 0.0 
            end
            #@show sourceₛ[ic],aₚ[ic],sourceᵤ[ic],sourceₚ[ic],t̂_dot_lf[ic,1:npc],ℾ,phi[ic]
            
        end

    end
   
    
    
end




function IC(cells,phi,ICvalue)
    for (ic,cell) in enumerate(eachrow(cells))
        phi[ic] = ICvalue
    end
    return nothing
end

function BC(boundary_edges,phi_bface,phi_node,BCarray)
    for (ib,bface) in enumerate(eachrow(boundary_edges))
        #@show i,bface
        for bvalue in BCarray
           # @show bface[5], bvalue.ID
            if bface[5] == bvalue.ID
                if bvalue.type == 1 # Dirichlet
                    phi_bface[ib] = bvalue.value
                end
                if bvalue.type == 2 # Neumann
                    phi_bface[ib] = (phi_node[bface[1]] + phi_node[bface[2]])/2
                end
            end
        end
    end
    return nothing
end




function main()
    BCarray = Array{BC, 1}(undef, 4)

    BCarray[1] = BC(5,1,1.0)
    BCarray[2] = BC(6,1,0.0)
    BCarray[3] = BC(7,2,0.0)
    BCarray[4] = BC(8,1,1.0)

    Sourceᵤ = 0.0
    Sourceₚ = 0.0

    ICvalue = 0.0

    ℾ = 0.1
    model_params = [ℾ]

    model_params = 0.0

    filename = ("864691136812081779.ply")
    
    uniquenode(filename,2,BCarray,ICvalue,ℾ,Sourceᵤ,Sourceₚ)

end



@time main()