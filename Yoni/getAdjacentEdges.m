function ZVector = getAdjacentEdges(nodeNum, edge, numOfNodes)
    
    ZVector(1) = nodeNum;
    ZCounter = 2;
    for i = 1:size(edge,1)
        if edge(i,1) == nodeNum
            ZVector(ZCounter) = i + numOfNodes; % Because after this, we're using this number to map to globalA.
            ZCounter = ZCounter + 1;
        end
        if edge(i,2) == nodeNum
            ZVector(ZCounter) = i + numOfNodes;
            ZCounter = ZCounter + 1;
        end
    end
    
end