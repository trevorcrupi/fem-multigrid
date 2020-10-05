function error = getL2Error(c,numOfTriangles,p,t)

    error          = 0;
    u              = @(x,y) (-cos(pi.*y)/pi);

    
    for i = 1:numOfTriangles      % Calculate (complete L2-error)^2
        columnVector = t(1:3, i); % Vector with points of the triangle
        localPhiCoeffs = getLocalPhiCoeffs(p,columnVector);
    
    
        phi1 = localPhiCoeffs(:,1);
        phi2 = localPhiCoeffs(:,2);
        phi3 = localPhiCoeffs(:,3);

        localU_h = @(x,y)                                              ...
                   c(columnVector(1)).*[phi1(1).*x+phi1(2).*y+phi1(3)] ...
                 + c(columnVector(2)).*[phi2(1).*x+phi2(2).*y+phi2(3)] ...
                 + c(columnVector(3)).*[phi3(1).*x+phi3(2).*y+phi3(3)];
    
             
        [X,Y,Wx,Wy]         = triquad(8, [p(1,columnVector(1)) p(2,columnVector(1)); p(1,columnVector(2)) p(2,columnVector(2)); p(1,columnVector(3)) p(2,columnVector(3))]); %triquad(8,[solvePhi(1,1) solvePhi(1,2); solvePhi(2,1) solvePhi(2,2); solvePhi(3,1) solvePhi(3,2)]);
        localL2ErrorSquared = Wx' * (u(X,Y)-localU_h(X,Y)).^2 * Wy;
        error               = error + localL2ErrorSquared;
     
    end
    error = sqrt(error);
end