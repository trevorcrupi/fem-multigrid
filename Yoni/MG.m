function solution = MG(inputUVector,globalB, meshNum, storingA, storingEdge, storeNodeNums,storeHeights)
% For "Test MG": add "MGErrorConvergenceRate,numOfMGIterations" as
%                more function outputs.


%     % Test MG
%     normS_iVector       = zeros(1,1);
%     errors3             = zeros(1,1);
%     normS_0             = getHrCurlNormforProblem3(inputUVector,storingA{meshNum});
%     checkForToleranceMG = 1;
%     numOfMGIterations   = 0;
    
% % % %     solution            = inputUVector;

%     while checkForToleranceMG >= (10^(-15))             % For "Test MG".
%         numOfMGIterations    = numOfMGIterations +1;

        if meshNum == 1 %Base case - ends recursion.
            solution   = storingA{meshNum}\globalB;
        
            
        else
            
            % Step 1
            V1 = GSFunctionForMG(storingA{meshNum}, globalB, storeHeights{meshNum}, storeNodeNums{meshNum}, inputUVector, meshNum);


            % Step 2
            zeroVector = zeros(storeHeights{meshNum-1},1);
            P          = prolong(meshNum-1);
            newGlobalB = P'*(globalB - storingA{meshNum} * V1);
            V2         = V1 + P * MG( zeroVector, newGlobalB, meshNum-1, storingA, storingEdge, storeNodeNums, storeHeights);


            % Step 3
            V3       = GSFunctionTranspose(storingA{meshNum}, globalB, storeHeights{meshNum}, storeNodeNums{meshNum}, V2, meshNum);
            solution = V3;
        end
    
        
%         % Test MG
%         if(numOfMGIterations == 1)
%            newval1                          = getHrCurlNormforProblem3(solution,storingA{meshNum});
%            normS_iVector(numOfMGIterations) = newval1;
%            errors3(numOfMGIterations)       = normS_iVector(numOfMGIterations) / normS_0;
%         end
%         if(numOfMGIterations > 1)
%             newval2       = getHrCurlNormforProblem3(solution,storingA{meshNum});
%             normS_iVector = [normS_iVector;newval2];
%             newval3       = normS_iVector(numOfMGIterations) / normS_iVector(numOfMGIterations-1);
%             errors3       = [errors3;newval3];
%         end
%         checkForToleranceMG = normS_iVector(numOfMGIterations) / normS_0;
        
%     end
    


%     % Test MG
%     %To check if this is working.
%     MGErrorConvergenceRate     = 0;
%     for l = 1:numOfMGIterations
%         MGErrorConvergenceRate = MGErrorConvergenceRate + errors3(l);
%     end
%     MGErrorConvergenceRate     = MGErrorConvergenceRate/numOfMGIterations;
    
end