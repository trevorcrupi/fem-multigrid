function [solution,MGErrorConvergenceRate,numOfMGIterations] = MGMain(inputUVector,globalB, meshNum, storingA, storingEdge, storeNodeNums,storeHeights)
% For "Test MG": add "MGErrorConvergenceRate,numOfMGIterations" as
%                more function outputs.


    % Test MG
    normS_iVector       = zeros(1,1);
    errors3             = zeros(1,1);
    normS_0             = getHrCurlNormforProblem3(inputUVector,storingA{meshNum});
    checkForToleranceMG = 1;
    numOfMGIterations   = 0;
    
    solution_i          = inputUVector;

    while checkForToleranceMG >= (10^(-15))             % This condition for "Test MG".
        numOfMGIterations    = numOfMGIterations +1;
        
        solution_i = MG( solution_i, globalB, meshNum, storingA, storingEdge, storeNodeNums,storeHeights);
    
        
        % Test MG
        if(numOfMGIterations == 1)
           newval1                          = getHrCurlNormforProblem3(solution_i,storingA{meshNum});
           normS_iVector(numOfMGIterations) = newval1;
           errors3(numOfMGIterations)       = normS_iVector(numOfMGIterations) / normS_0;
        end
        if(numOfMGIterations > 1)
            newval2       = getHrCurlNormforProblem3(solution_i,storingA{meshNum});
            normS_iVector = [normS_iVector;newval2];
            newval3       = normS_iVector(numOfMGIterations) / normS_iVector(numOfMGIterations-1);
            errors3       = [errors3;newval3];
        end
        checkForToleranceMG = normS_iVector(numOfMGIterations) / normS_0;
    end
    solution = solution_i;
    


    % Test MG
    MGErrorConvergenceRate     = 0;
    for l = 1:numOfMGIterations
        MGErrorConvergenceRate = MGErrorConvergenceRate + errors3(l);
    end
    MGErrorConvergenceRate     = MGErrorConvergenceRate/numOfMGIterations;
end