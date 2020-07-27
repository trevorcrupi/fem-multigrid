clear all

iterations              = 3; %Enough PETs and new_Eles to do 8.
meshLevel               = zeros(iterations,1);
sizeOfDim               = zeros(iterations,1);   %Test GS1
errorConvergenceRates   = zeros(iterations,1);   %Test GS1
numOfGSIterations       = zeros(iterations,1);   %Test GS1
k                       = 1; % Pick any non-zero integer, positive or negative. 1, -1, 2, -1; 10-(-10) For now.

for n = 1:iterations
    
    mystr             = ['PETForYoni/PETForYoni' num2str(n) '.mat'];
    load(mystr);
    edge              = getEdgeMatrix(p,t);
    storingEdge{n}    = edge;
    meshLevel(n)      = n;
    numOfTriangles(n) = size(t,2);

    %Test GS1
    [c,errorConvergenceRate,numOfGSIterations1,numOfNodesAndEdges] = solveApproximationForTestGS(p,e,t,numOfTriangles(n),k,edge,n);
    errorConvergenceRates(n)                                       = errorConvergenceRate;
    numOfGSIterations(n)                                           = numOfGSIterations1;
    sizeOfDim(n)                                                   = numOfNodesAndEdges;
    
end

% Test GS1 . . .
MeshLevel            = zeros(iterations,1);
SizeOfDim            = zeros(iterations,1);
NumOfGSIterations    = zeros(iterations,1);
ErrorConvergenceRate = zeros(iterations,1);

for j = 1:iterations
    MeshLevel(j)            = meshLevel(j);
    SizeOfDim(j)            = sizeOfDim(j);
    NumOfGSIterations(j)    = numOfGSIterations(j);
    ErrorConvergenceRate(j) = errorConvergenceRates(j);
end

table(MeshLevel,SizeOfDim,NumOfGSIterations,ErrorConvergenceRate)