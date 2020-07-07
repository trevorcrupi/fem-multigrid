function localPhiCoeffs = getLocalPhiCoeffs(p,columnVector)

    solvePhi = [p(1,columnVector(1)) p(2,columnVector(1)) 1
                p(1,columnVector(2)) p(2,columnVector(2)) 1
                p(1,columnVector(3)) p(2,columnVector(3)) 1];

    i1 = [1;0;0];
    i2 = [0;1;0];
    i3 = [0;0;1];

    % Matrix with columns of Phi Coefficients
    localPhiCoeffs = [ linsolve(solvePhi, i1) linsolve(solvePhi, i2) linsolve(solvePhi, i3) ];

end