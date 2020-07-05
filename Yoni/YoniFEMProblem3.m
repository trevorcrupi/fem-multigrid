clear all
load PETForYoni/PETForYoni1.mat
load newEle/new_ele1.mat

edge = getEdgeMatrix(p,t);

numOfTriangles = size(t,2);
 for i = 1:numOfTriangles

    columnVector = t(1:3, i); % Vector with points of the triangle
    rowVector = new_ele(i, 1:3); % Vector with points of the triangle
    
    localPhiCoeffs = getLocalPhiCoeffs(p,columnVector);
    localPsiCoeffs = getLocalPsiCoeffs(p,rowVector,edge);
    
 end