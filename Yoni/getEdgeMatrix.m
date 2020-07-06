function edge = getEdgeMatrix(p,t)

    [~,N_node]=size(p); %N_node: number of nodes
    node=p';  %node: the N_node by 2 matrix that has the x,y-coordinates of each node.
    [~,N_ele]=size(t); %N_ele is the number of triangles.
    ele=t(1:3,1:N_ele);
    ele=ele'; %ele: the N_ele by 3 matrix that has the three vertex numbers for each triangle.
    TR=triangulation(ele,node);
    edge=edges(TR); %edge: the N_edge by 2 matrix that has the two vertex numbers for making each edge.
    [N_edge,E2]=size(edge);  %N_edge is the number of edges.

end