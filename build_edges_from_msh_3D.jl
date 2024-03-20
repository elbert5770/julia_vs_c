function edge_comparison(edge_to_cell_sort,edges,boundary_edges,cells,node_list,npc,count1,count2,skip)

    for (i,ftc) in enumerate(eachrow(edge_to_cell_sort))
        if skip == 1
            node_list[1,1:5] = ftc
            skip = 0
            continue
        end
        
        node_list[2,1:5] = ftc
        if node_list[1,1:2] == node_list[2,1:2]
            count1 += 1

            # edges = [node 1 of both cells, node 2 of both cells, cell 1,cell 2, zone,boundary edge?,edge# of edge in cell 1,edge# of edge in cell 2,boundary edge #,node order cell 1 reversed?,node order cell 2 reversed?]
            # node_list = [node 1 of both cells, node 2 of both cells, cell 1, edge# of edge in cell 1 (1 to npc),was order of nodes reversed?
            #              node 1 of both cells, node 2 of both cells, cell 2, edge# of edge in cell 1 (1 to npc),was order of nodes reversed?]
            # order of nodes reveresed? 0 = false, 1 = true; order needs to be reversed for sorting and matching
            edges[count1,1:11] = [node_list[1,1],node_list[1,2],node_list[1,3],node_list[2,3],1,0,node_list[1,4],node_list[2,4],0,node_list[1,5],node_list[2,5]]
            cells[node_list[1,3],npc+node_list[1,4]] = count1
            cells[node_list[2,3],npc+node_list[2,4]] = count1
            cells[node_list[1,3],npc*2+node_list[1,4]] = node_list[2,3]
            cells[node_list[2,3],npc*2+node_list[2,4]] = node_list[1,3]
            cells[node_list[1,3],npc*3+1] = 1
            cells[node_list[2,3],npc*3+1] = 2
            skip = 1
        else
            count1 += 1
            count2 += 1
      
            edges[count1,1:11] = [node_list[1,1],node_list[1,2],node_list[1,3],0,5,1,node_list[1,4],0,count2,node_list[1,5],0]
            boundary_edges[count2,1:11] = edges[count1,1:11]
            cells[node_list[1,3],npc+node_list[1,4]] = count1
            cells[node_list[1,3],npc*2+node_list[1,4]] = 0
            

            node_list[1,1:5] = node_list[2,1:5]
            skip = 0
        end
    end
    if skip == 0
        count1 += 1
        count2 += 1
        edges[count1,1:11] = [node_list[1,1],node_list[1,2],node_list[1,3],0,5,1,node_list[1,4],0,count2,node_list[1,5],0]
        boundary_edges[count2,1:11] = edges[count1,1:11]
        cells[node_list[1,3],npc+node_list[1,4]] = count1
        cells[node_list[1,3],npc*2+node_list[1,4]] = 0
    end
    return count1,count2
end



function link_cells_edges(edges,cells,num_edges,npc)

    cell_temp = Array{Int64,2}(undef, num_edges*2,5)

    
    for (ie,edge) in enumerate(eachrow(edges))
        node1 = edge[1]
        node2 = edge[2]
        cell1 = edge[3]
        cell2 = edge[4]
        zone = edge[5]

        if node1 < node2
            cell_temp[(ie-1)*2+1,1:3] = [cell1,node1,node2]
            cell_temp[(ie-1)*2+2,1:3] = [cell2,node1,node2]
        else
            cell_temp[(ie-1)*2+1,1:3] = [cell1,node2,node1]
            cell_temp[(ie-1)*2+2,1:3] = [cell2,node2,node1]       
        end   
        cell_temp[(ie-1)*2+1,4] = zone
        cell_temp[(ie-1)*2+2,4] = zone
        cell_temp[(ie-1)*2+1,5] = ie
        cell_temp[(ie-1)*2+2,5] = ie

        # Take advantage of the loop to set default boundary edge flag to zero, corresponding to not a boundary edge
        # edges[ie,6] = 0
    end
    
    cell_temp = sortslices(cell_temp,dims=1,by=x->(x[1],x[2],x[3]))
   # @show nodes[2213:2220,:]
    @show cell_temp[1:100,:]
    @show cell_temp[end-10:end,:]
    
    count5 = 0
    for cellt in eachrow(cell_temp)

        # Each cell has three edges.  Each cell's three edges are stored in
        # three rows of cell_temp.  These three rows are combined into one
        # row of 'cells', with one row per cell.  The order of the nodes
        # and edges is in ascending order of nodes, so the direction
        # around the triangle that the nodes and edges are ordered is 
        # not fixed.
        # cell has the following fields:
        # [zone,node1,node2,node3,edge12,edge23,edge31]
        # edge1 is between node1 & node2
        # edge2 is between node2 & node3
        # edge3 is between node3 & node1

        cell1 = cellt[1]
        node1 = cellt[2]
        node2 = cellt[3]
        zone = cellt[4]
        edge12 = cellt[5]
        if cell1 == 0
            # edges[edge12,6] = 1
            continue
        end
                
        if npc == 4
            if mod(count5,npc) == 0
                cells[cell1,1] = node1
                cells[cell1,2] = node2
                cells[cell1,5] = edge12
            else
                if mod(count5,npc) == 2
                    if cells[cell1,2] == node1                   
                        cells[cell1,3] = node2               
                        cells[cell1,6] = edge12
                    else
                        cells[cell1,3] = node1               
                        cells[cell1,6] = edge12
                    end
                else
                    if mod(count5,npc) == 1
                        cells[cell1,8] = edge12
                        cells[cell1,4] = node2
                    else
                        cells[cell1,7] = edge12
                        
                    end

                end
            end
            count5 += 1
        end

        if npc == 3
            if mod(count5,npc) == 0
                cells[cell1,1] = node1
                cells[cell1,2] = node2
                cells[cell1,4] = edge12
            else
                if mod(count5,npc) == 2
                    if cells[cell1,2] == node1                   
                        cells[cell1,3] = node2               
                        cells[cell1,5] = edge12
                    else
                        cells[cell1,3] = node1               
                        cells[cell1,5] = edge12
                    end
                else
                    if mod(count5,npc) == 1
                        cells[cell1,6] = edge12
        
                    end

                end
            end
            count5 += 1
        end


    end
    for (ic,cell) in enumerate(eachrow(cells))

        for j in range(1,npc)
            if ic == edges[cell[j+npc],3]
                cells[ic,j+npc*2] = edges[cell[j+npc],4]
            else
                cells[ic,j+npc*2] = edges[cell[j+npc],3]
            end
        end
    end
    @show cells[1:10,:]
    @show cells[end-10:end,:]

    
    


end


function fix_edges()
    num_bfaces = sum(edges[:,6])
    
    area_face = Vector{Float64}(undef,num_edges)
    
    bfaces = Array{Int64,2}(undef,num_bfaces,5)
   
    area_bface = Vector{Float64}(undef,num_bfaces)
    
    count2 = 0

    for (ie,edge) in enumerate(eachrow(edges))
        node1 = edge[1]
        node2 = edge[2]
        cell1 = edge[3]
        cell2 = edge[4]
        boundary_face = edge[6]


        area_face[ie] = norm(nodes[node2,:]-nodes[node1,:])
        
        if boundary_face == 0 # not a boundary face

            # Record in edges which of the cell's three faces is the edge
            # Will be a number between 1 and 3.
            for j in range(1,npc)
                if cells[cell1,j+npc] == ie
                    edges[ie,7] = j
                end
                if cells[cell2,j+npc] == ie 
                    edges[ie,8] = j
                end
            end
    
           
            
        else
            count2 += 1
            bfaces[count2,1:5] = edge[1:5]
           
            for j in range(1,npc)
                if cells[cell1,j+npc] == ie 
                    bfaces[count2,4] = j
                    edges[ie,7] = j
                end
            end
            edges[ie,8] = 0 #  there is no second cell to have a face number
            edges[ie,9] = count2

            area_bface[count2] = area_face[ie]
        end

    end
end