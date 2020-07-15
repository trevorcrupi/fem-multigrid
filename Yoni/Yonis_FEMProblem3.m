clear all

iterations     = 6;
meshLevel      = zeros(iterations,1);
numOfTriangles = zeros(iterations,1);
errorVector    = zeros(iterations,1);
rateOfErrors   = zeros(iterations-1,1);
k              = 1; % Pick any non-zero integer, positive or negative. 1, -1, 2, -1; 10-(-10) For now.

for n = 1:iterations
    
    mystr = ['PETForYoni/PETForYoni' num2str(n) '.mat'];
    load(mystr);
    edge           = getEdgeMatrix(p,t);
    meshLevel(n)   = n;
    
    numOfTriangles(n) = size(t,2);
    c                 = solveApproximationForProb3(p,e,t,numOfTriangles(n),k,edge,n);
    
    % Before Multigrid
    errorVector(n)    = getL2ErrorforProblem3(c,numOfTriangles(n),p,t,edge,k,n);
    if(n > 1)
        rateOfErrors(n-1) = log2( errorVector(n-1) / errorVector(n) );
    end

    %With GS1
%     errorVector(n)    = getL2ErrorforProblem3(c,numOfTriangles(n),p,t,edge,k,n);
%     if(n > 1)
%         rateOfErrors(n-1) = ( errorVector(n) / errorVector(n-1) );
%     end
end

MeshLevel      = zeros(iterations-1,1);
NumOfTriangles = zeros(iterations-1,1);
ErrorVector    = zeros(iterations-1,1);
RateOfErrors   = zeros(iterations-1,1);

for j=1:iterations-1
    MeshLevel(j)      = meshLevel(j);
    NumOfTriangles(j) = numOfTriangles(j);
    ErrorVector(j)    = errorVector(j);
    RateOfErrors(j)   = rateOfErrors(j);
end

table(MeshLevel,NumOfTriangles,ErrorVector,RateOfErrors)