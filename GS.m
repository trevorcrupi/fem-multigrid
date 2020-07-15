

function solution = GS(globalA, globalB, edge, N)
    M = size(edge, 1);
    currentU = ones(N+M, 1);
    OGError = sqrt(currentU' * globalA * currentU);
    
    for i=1:N
        Z_i = getAdjacentEdges(i, N, edge);
        ZDim = length(Z_i);
        A_j = zeros(ZDim, ZDim);
        fminusAU = globalB - globalA*currentU;
        
        projection = zeros(ZDim, 1);
        for j=1:ZDim
            for k=1:ZDim
                A_j(j, k) = globalA(Z_i(j), Z_i(k));
            end 
            
            projection(j) = fminusAU(Z_i(j));
        end 
        
        inverseATimesProjectionFullVector = zeros(N+M, 1);
        inverseATimesProjection = inv(A_j)*projection;
        for o=1:ZDim
            inverseATimesProjectionFullVector(Z_i(o)) = inverseATimesProjection(o);
        end 
        
        currentU = currentU + inverseATimesProjectionFullVector;
        newError = sqrt(currentU'*globalA*currentU);
        
        newError/OGError
    end
    
    solution = currentU;
end