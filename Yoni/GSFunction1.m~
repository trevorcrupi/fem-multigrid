function [solution,errorConvergenceRate,numOfGSIterations] = GSFunction1(globalA, globalB, height, numOfNodes, edge, inputUVector)

    U          = inputUVector;
    normU_iVector = zeros(1,1);
    errors2  = zeros(1,1);
    norm_0         = getHrCurlErrorforProblem3(inputUVector,globalA);
    %checkForTolerance = zeros(iterations,1);
    checkForTolerance = 1;
    numOfGSIterations = 0;
    U_i = zeros(size(inputUVector,1));
    
    
    while checkForTolerance >= (10^(-15))
        
        numOfGSIterations = numOfGSIterations +1;
        for i = 1:numOfNodes

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



        end
        U_i = U_i + U;

        
        if(numOfGSIterations == 1)
           normU_iVector(numOfGSIterations) = getHrCurlErrorforProblem3(U_i,globalA);
           errors2(numOfGSIterations) = normU_iVector(numOfGSIterations) / norm_0;
        end
        if(numOfGSIterations > 1)
            normU_iVector = [normU_iVector;getHrCurlErrorforProblem3(U_i,globalA)];
            newval = normU_iVector(numOfGSIterations) / normU_iVector(numOfGSIterations-1);
            errors2 = [errors2;newval];
        end
        
        checkForTolerance = normU_iVector(numOfGSIterations) / norm_0;
    end
    
    solution = U_i;
    numOfGSIterations;
    
    
    %To check if this is working.
    errorConvergenceRate = 0;
    for l = 1:numOfGSIterations
        errorConvergenceRate = errorConvergenceRate + errors2(l);
    end
    errorConvergenceRate = errorConvergenceRate/numOfGSIterations;
    
end