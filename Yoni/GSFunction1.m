function [solution,errorConvergenceRate] = GSFunction1(globalA, globalB, height, numOfNodes, edge, inputUVector)

    iterations = numOfNodes;
    U          = inputUVector;
    normU_iVector = zeros(iterations,1);
    errors2  = zeros(iterations,1);
    norm_0         = getHrCurlErrorforProblem3(inputUVector,globalA);
    %checkForTolerance = zeros(iterations,1);
    checkForTolerance = 1;
    i = 0;
    
    
    while (checkForTolerance >= (10^(-15))) && (i <= numOfNodes)
    %for i = 1:iterations
        i = i +1;
        
        Z_i      = getAdjacentEdges(i, edge, numOfNodes);
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
        
        
        
        
        
        
        normU_iVector(i) = getHrCurlErrorforProblem3(U,globalA);
        
        if(i == 1)
           errors2(i) = ( normU_iVector(i) / norm_0 );
        end
        if(i > 1)
            errors2(i) = ( normU_iVector(i) / normU_iVector(i-1) );
        end
        
        checkForTolerance = normU_iVector(i) / norm_0;
        
    %end
    end
    
    
    solution = U;
    
    
    
    errorConvergenceRate = 0;
    for l = 1:i
        errorConvergenceRate = errorConvergenceRate + errors2(l);
    end
    errorConvergenceRate = errorConvergenceRate/i;
    i
    
end