function edges = getAdjacentEdges(node, N, edge) 
    j = 2;
    edges(1) = node;
    for i=1:size(edge, 1)
        if edge(i, 1) == node
            edges(j) = i+N;
            j = j+1;
        end 
        if edge(i, 2) == node
            edges(j) = i+N;
            j = j+1;
        end
        
    end
end