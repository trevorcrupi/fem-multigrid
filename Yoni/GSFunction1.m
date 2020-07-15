function solution = GSFunction1(globalA, globalB, height, numOfNodes, edge, inputUVector)

    U        = inputUVector;
    
    for i = 1:numOfNodes
        
        Z_i      = getAdjacentEdges(i, edge, numOfNodes)  ;
        ZDim     = length(Z_i);
        A_j      = zeros(ZDim, ZDim);
        fMinusAU = globalB - globalA * U ;
        
        
        % Fill in A_j
        for z = 1:ZDim
            for y = 1:ZDim
                A_j(z,y) = globalA(Z_i(z), Z_i(y));
            end
        end
        
        
        fMinusAUAfterQ = zeros(ZDim,1);
        for j = 1:ZDim
            fMinusAUAfterQ(j) = fMinusAU(Z_i(j));
        end
        
        
        invATimesProjection           = inv(A_j) * fMinusAUAfterQ;
        invATimesProjectionFullVector = zeros(height,1);
        for o = 1:ZDim
            invATimesProjectionFullVector(Z_i(o)) = invATimesProjection(o);
        end
        
        
        U = U + invATimesProjectionFullVector;
    end
    
    solution = U;
    
end