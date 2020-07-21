load PETForYoni/PETForYoni1.mat
load newEle/new_ele1.mat

numOfNodes = size(p,2);
edge = getEdgeMatrix(p,t);

for i = 1:numOfNodes
    getAdjacentEdges(i,edge,numOfNodes)
end