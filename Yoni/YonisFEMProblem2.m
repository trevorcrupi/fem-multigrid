load PETForYoni/PETForYoni1.mat

numOfTriangles = size(t,2);
for i = 1:numOfTriangles
    
    columnVector = t(1:3, i); % Vector with points of the triangle
    
    % Points of triangle
    p1 = columnVector(1);
    p2 = columnVector(2);
    p3 = columnVector(3);
        
        
    solvePsi = [ p(1,p1).*p(2,p2)-p(2,p1).*p(1,p2)   p(1,p2)-p(1,p1)   p(2,p2)-p(2,p1)            %[ x1y2-y1x2  x2-x1  y2-y1
                 p(1,p2).*p(2,p3)-p(2,p2).*p(1,p3)   p(1,p3)-p(1,p2)   p(2,p3)-p(2,p2)            %  x1y2-y1x2  x2-x1  y2-y1
                 p(1,p1).*p(2,p3)-p(2,p1).*p(1,p3)   p(1,p3)-p(1,p1)   p(2,p3)-p(2,p1) ];         %  x1y2-y1x2  x2-x1  y2-y1 ]

    i1 = [1;0;0];
    i2 = [0;1;0];
    i3 = [0;0;1];

    % Matrix with columns of Phi Coefficients
    localPsiCoeffs = [ linsolve(solvePsi, i1) linsolve(solvePsi, i2) linsolve(solvePsi, i3) ];
    
    tangent = [ p(1,p3)-p(1,p2); p(2,p3)-p(2,p2) ].*sqrt((p(1,p3)-p(1,p2)).^2+(p(2,p3)-p(2,p2)).^2);
end