function error = getL2ErrorSquaredforProblem3(c,numOfTriangles,p,t)

    error          = 0;
    
    U_r            = @(r,z) z - (1/k).*( (r.^3)/3 - (r.^2)/2 );
    U_theta        = @(r,z) -k.*z + (r.^3)/3 - (r.^2)/2;
    U_z            = @(r,z) r;
    
    
    for i = 1:numOfTriangles      % Calculate (complete L2-error)^2
        
        columnVector     = t(1:3, i); % Vector with points of the triangle
        rowVector        = new_ele(i, 1:3); % Vector with edges of the triangle
        localPhiCoeffs   = getLocalPhiCoeffs(p,columnVector);
        localPsiCoeffs   = getLocalPsiCoeffs(p,rowVector,edge);
        updatedRowVector = rowVector + numOfNodes; % Now numbered for n-nodes + N-edges -- bumps up numbering of edges, so node numbers are first.
        cAndRVector      = [columnVector; updatedRowVector'];
    
        phi1 = localPhiCoeffs(:,1);
        phi2 = localPhiCoeffs(:,2);
        phi3 = localPhiCoeffs(:,3);
        psi1 = localPsiCoeffs(:,1);
        psi2 = localPsiCoeffs(:,2);
        psi3 = localPsiCoeffs(:,3);
        

        localU_h = @(x,y)                                              ...
                   c(cAndRVector(1)).*[phi1(1).*x+phi1(2).*y+phi1(3)] ...
                 + c(cAndRVector(2)).*[phi2(1).*x+phi2(2).*y+phi2(3)] ...
                 + c(cAndRVector(3)).*[phi3(1).*x+phi3(2).*y+phi3(3)];
    
             
        [X,Y,Wx,Wy]         = triquad(8, [p(1,columnVector(1)) p(2,columnVector(1)); p(1,columnVector(2)) p(2,columnVector(2)); p(1,columnVector(3)) p(2,columnVector(3))]); %triquad(8,[solvePhi(1,1) solvePhi(1,2); solvePhi(2,1) solvePhi(2,2); solvePhi(3,1) solvePhi(3,2)]);
        localL2ErrorSquared = Wx' * ( r.* [u(X,Y) - localU_h(X,Y)] ).^2 * Wy;
        
        % U_r - U_theta - U_z -
        
        
        error               = error + localL2ErrorSquared;
     
    end
    error = sqrt(error);
end