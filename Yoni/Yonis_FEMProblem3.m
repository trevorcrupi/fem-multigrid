clear all

iterations              = 8; %Enough PETs and new_Eles to do 8.
meshLevel               = zeros(iterations,1); 
numOfTriangles          = zeros(iterations,1);   %Before Multigrid
errorVector             = zeros(iterations,1);   %Before Multigrid
rateOfErrors            = zeros(iterations-1,1); %Before Multigrid
sizeOfDim               = zeros(iterations,1);   %Test GS1
errorConvergenceRates   = zeros(iterations,1);   %Test GS1
numOfGSIterations       = zeros(iterations,1);   %Test GS1
storingA                = cell(1,8); % Put in globalA.    %With Multigrid (Test MG)
storingEdge             = cell(1,8); % Put in edge.       %With Multigrid (Test MG)
storeNodeNums           = cell(1,8); % Put in numOfNodes. %With Multigrid (Test MG)
storeHeights            = cell(1,8); % Put in height.     %With Multigrid (Test MG)
errorConvergenceRatesMG = zeros(iterations,1);            %With Multigrid (Test MG)
numOfMGIterations       = zeros(iterations,1);            %With Multigrid (Test MG)
sizeOfDimForMG          = zeros(iterations,1);            %With Multigrid (Test MG)
k                       = 1; % Pick any non-zero integer, positive or negative. 1, -1, 2, -1; 10-(-10) For now.

for n = 1:iterations
    
    mystr             = ['PETForYoni/PETForYoni' num2str(n) '.mat'];
    load(mystr);
    edge              = getEdgeMatrix(p,t);
    storingEdge{n}    = edge;
    meshLevel(n)      = n;
    numOfTriangles(n) = size(t,2);
    
    %Before Multigrid
%     c                     = solveApproximationForProb3(p,e,t,numOfTriangles(n),k,edge,n);
%     errorVector(n)        = getL2ErrorforProblem3(c,numOfTriangles(n),p,t,edge,k,n);
%     if(n > 1)
%         rateOfErrors(n-1) = log2( errorVector(n-1) / errorVector(n) );
%     end

    %Test GS1
    [c,errorConvergenceRate,numOfGSIterations1,numOfNodesAndEdges] = solveApproximationForProb3(p,e,t,numOfTriangles(n),k,edge,n);
    errorConvergenceRates(n)                                       = errorConvergenceRate;
    numOfGSIterations(n)                                           = numOfGSIterations1;
    sizeOfDim(n)                                                   = numOfNodesAndEdges;
    
    %With Multigrid (Test MG)
%     [c,globalA,height,numOfNodes,MGErrorConvergenceRate,numOfMGIterations1] = solveApproximationForProb3(p,e,t,numOfTriangles(n),k,edge,n,storingA,storingEdge,storeNodeNums,storeHeights);
% %     MGErrorConvergenceRate,numOfMGIterations1
%     storingA{n}                = globalA;
%     storeNodeNums{n}           = numOfNodes;
%     storeHeights{n}            = height;
%     sizeOfDimForMG(n)          = height;
%     errorConvergenceRatesMG(n) = MGErrorConvergenceRate;
%     numOfMGIterations(n)       = numOfMGIterations1;
    
end

% Before Multigrid . . .
% MeshLevel      = zeros(iterations-1,1);
% NumOfTriangles = zeros(iterations-1,1);
% ErrorVector    = zeros(iterations-1,1);
% RateOfErrors   = zeros(iterations-1,1);
% 
% for j=1:iterations-1
%     MeshLevel(j)      = meshLevel(j);
%     NumOfTriangles(j) = numOfTriangles(j);
%     ErrorVector(j)    = errorVector(j);
%     RateOfErrors(j)   = rateOfErrors(j);
% end
% 
% table(MeshLevel,NumOfTriangles,ErrorVector,RateOfErrors)

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

% With Multigrid (Test MG) . . .
% MeshLevel                 = zeros(iterations,1);
% SizeOfDim                 = zeros(iterations,1);
% NumOfMGIterations         = zeros(iterations,1);
% ErrorConvergenceRateForMG = zeros(iterations,1);
% 
% for j = 1:iterations
%     MeshLevel(j)                 = meshLevel(j);
%     SizeOfDimForMG(j)            = sizeOfDimForMG(j);
%     NumOfMGIterations(j)         = numOfMGIterations(j);
%     ErrorConvergenceRateForMG(j) = errorConvergenceRatesMG(j);
% end
% 
% table(MeshLevel,SizeOfDimForMG,NumOfMGIterations,ErrorConvergenceRateForMG)