function solution = GSFunctionTranspose(globalA, globalB, height, numOfNodes, inputUVector, meshNum)
% For "With Multigrid": Runs for 1 iteration.


    U_i               = inputUVector;
    mystr             = ['ZMatrix/ZMatrix' num2str(meshNum) '.mat'];
    load(mystr);
    
    
    for i = numOfNodes:1
        
        Z_i                                = ZMatrix(i,:);
        Z_i                                = Z_i(Z_i > 0);
        diff                               = globalB - globalA * U_i;
        diff                               = diff(Z_i);
        A_j                                = globalA(Z_i,Z_i);
        invATimesProjection                = A_j\diff;
        invATimesProjectionFullVector      = zeros(height,1);
        invATimesProjectionFullVector(Z_i) = invATimesProjection;
        U_i                                = U_i + invATimesProjectionFullVector;
    
    end
    solution = U_i;
end