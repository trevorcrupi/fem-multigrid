function error = getL2ErrorSquaredforProblem3(c,numOfTriangles,p,t,edge,k,iterationNum)

    mystr = ['newEle/new_ele' num2str(iterationNum) '.mat'];
    load(mystr);
    
    numOfNodes     = size(p,2);
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
        
        phiLine1 = @(basis,r,z) [ - (basis(3)/k) - (basis(1)/k).*r - (basis(2)/k).*z ];
        phiLine2 = @(basis,r,z) [ basis(1).*r + basis(2).*z + basis(3) ];
        phiLine3 = @() 0;
        psiLine1 = @(basis,r,z) [ (basis(2)/k).*r - (basis(1)/k).*r.*z ];
        psiLine2 = @() 0;
        psiLine3 = @(basis,r) [ (basis(3)/k).*r + (basis(1)/k).*r.^2 ];
        

        localU_hLine1 = @(r,z) ...
                        c(cAndRVector(1))*phiLine1(phi1,r,z) + c(cAndRVector(2))*phiLine1(phi2,r,z) + c(cAndRVector(3))*phiLine1(phi3,r,z) ...
                      + c(cAndRVector(4))*psiLine1(psi1,r,z) + c(cAndRVector(5))*psiLine1(psi2,r,z) + c(cAndRVector(6))*psiLine1(psi3,r,z);
        localU_hLine2 = @(r,z) ...
                        c(cAndRVector(1))*phiLine2(phi1,r,z) + c(cAndRVector(2))*phiLine2(phi2,r,z) + c(cAndRVector(3))*phiLine2(phi3,r,z) ...
                      + c(cAndRVector(4))*psiLine2()         + c(cAndRVector(5))*psiLine2()         + c(cAndRVector(6))*psiLine2();
        localU_hLine3 = @(r,z) ...
                        c(cAndRVector(1))*phiLine3()       + c(cAndRVector(2))*phiLine3()       + c(cAndRVector(3))*phiLine3()       ...
                      + c(cAndRVector(4))*psiLine3(psi1,r) + c(cAndRVector(5))*psiLine3(psi2,r) + c(cAndRVector(6))*psiLine3(psi3,r);
    

        [X,Y,Wx,Wy]         = triquad(8, [p(1,columnVector(1)) p(2,columnVector(1)); p(1,columnVector(2)) p(2,columnVector(2)); p(1,columnVector(3)) p(2,columnVector(3))]);
        
        errorLine1 = U_r(X,Y)     - localU_hLine1(X,Y);
        errorLine2 = U_theta(X,Y) - localU_hLine2(X,Y);
        errorLine3 = U_z(X,Y)     - localU_hLine3(X,Y);
        
        localL2ErrorSquared = Wx' * (X .* [errorLine1.^2 + errorLine2.^2 + errorLine3.^2]) * Wy;
        error               = error + localL2ErrorSquared;
     
    end
    error = sqrt(error);
end