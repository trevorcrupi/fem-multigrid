function localPsiCoeffs = getLocalPsiCoeffs(p,rowVector,edge)
   
    edge1 = rowVector(1);
    edge2 = rowVector(2);
    edge3 = rowVector(3);
    
    % x1y2-y1x2  x2-x1  y2-y1
    solvePsi = [ p(1,edge(edge1,1)).*p(2,edge(edge1,2))-p(2,edge(edge1,1)).*p(1,edge(edge1,2))   p(1,edge(edge1,2))-p(1,edge(edge1,1))   p(2,edge(edge1,2))-p(2,edge(edge1,1))
                 p(1,edge(edge2,1)).*p(2,edge(edge2,2))-p(2,edge(edge2,1)).*p(1,edge(edge2,2))   p(1,edge(edge2,2))-p(1,edge(edge2,1))   p(2,edge(edge2,2))-p(2,edge(edge2,1))
                 p(1,edge(edge3,1)).*p(2,edge(edge3,2))-p(2,edge(edge3,1)).*p(1,edge(edge3,2))   p(1,edge(edge3,2))-p(1,edge(edge3,1))   p(2,edge(edge3,2))-p(2,edge(edge3,1)) ];
    
    
    i1 = [1;0;0];
    i2 = [0;1;0];
    i3 = [0;0;1];
 
    % Matrix with columns of Phi Coefficients
    localPsiCoeffs = [ linsolve(solvePsi, i1) linsolve(solvePsi, i2) linsolve(solvePsi, i3) ];
    
end