function integrandFunct = getIntegrand(basis_i, basis_j, basisI, basisJ, r, z, idNum)

    %k = ;
    curlRZKforPhi = @(basis) [   basis(2)
                           - basis(2)/k
                             basis(1)   ];
    curlRZKforPsi = @(basis,r,z) [ - basis(3) - basis(1).*r
                               - basis(3)/k - 3.*(basis(1)/k).*r
                                 basis(2) - basis(1).*z          ];

                             
                             
                             
    % *? or .*? for final r...
    
    if idNum = 1
        integrandFunct = [ ( curlRZKforPhi(basisI) * curlRZKforPhi(basisJ) + basis_i * basis_j) *r];
    end
    
    if idNum = 2
        integrandFunct = [ ( curlRZKforPhi(basisI) * curlRZKforPsi(basisJ) + basis_i * basis_j) *r];
    end
    
    if idNum = 3
        integrandFunct = [ ( curlRZKforPsi(basisI) * curlRZKforPhi(basisJ) + basis_i * basis_j) *r];
    end
    
    if idNum = 4
        integrandFunct = [ ( curlRZKforPsi(basisI) * curlRZKforPsi(basisJ) + basis_i * basis_j) *r];
    end
    
end