clear all
load PETForYoni/PETForYoni1.mat
load newEle/new_ele1.mat

edge = getEdgeMatrix(p,t);

% ourF = ;
curlRZK = @() [
    ];

psiVector = @(basis,r,z) [ (basis(2)/k).*r - (basis(1)/k).*r.*z
                            0
                           (basis(3)/k).*r + (basis(1)/k).*r.^2 ];
                       
phiVector = @(basis,r,z) [ - (basis(3)/k) - (basis(1)/k).*r - (basis(2)/k).*z
                           basis(1).*r + basis(2).*z + basis(3)
                                                               0];

                                                           
                                                           
numOfTriangles = size(t,2);
 for i = 1:numOfTriangles

    columnVector = t(1:3, i); % Vector with points of the triangle
    rowVector = new_ele(i, 1:3); % Vector with points of the triangle
    
    localPhiCoeffs = getLocalPhiCoeffs(p,columnVector);
    localPsiCoeffs = getLocalPsiCoeffs(p,rowVector,edge);
    
    
    
 end