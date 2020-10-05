clear all

iterations              = 6; %Enough PETs and new_Eles to do 8.
meshLevel               = zeros(iterations,1);
storingA                = cell(1,8); % Put in globalA.    %With Multigrid (Test MG)
storingEdge             = cell(1,8); % Put in edge.       %With Multigrid (Test MG)
storeNodeNums           = cell(1,8); % Put in numOfNodes. %With Multigrid (Test MG)
storeHeights            = cell(1,8); % Put in height.     %With Multigrid (Test MG)
errorConvergenceRatesMG = zeros(iterations,1);            %With Multigrid (Test MG)
numOfMGIterations       = zeros(iterations,1);            %With Multigrid (Test MG)
sizeOfDimForMG          = zeros(iterations,1);            %With Multigrid (Test MG)
k                       = 1; % Pick any non-zero integer, positive or negative. 1, -1, 2, -1; 10-(-10) For now.


%The new code . . .
load MyDomain.mat; %Load MyDomain file
g = decsg(gd, sf, ns);
model = createpde(1);
geometryFromEdges(model, g);
[p, e, t] = initmesh(g, 'hmax', inf);

for n = 1:iterations
    
    %The code that was replaced . . .
%     mystr             = ['PETForYoni/PETForYoni' num2str(n) '.mat'];
%     load(mystr);
    

    %The new code . . .
    if n > 1
        [p, e, t] = refinemesh(g, p, e, t, 'regular');
    end
    
    edge              = getEdgeMatrix(p,t);
    storingEdge{n}    = edge;
    meshLevel(n)      = n;
    numOfTriangles(n) = size(t,2);

% [c,globalA,height,numOfNodes] = solveApproximationForMultigrid(p,e,t,numOfTriangles(n),k,edge,n,storingA,storingEdge,storeNodeNums,storeHeights);
    %With Multigrid (Test MG)
    [c,globalA,height,numOfNodes,MGErrorConvergenceRate,numOfMGIterations1] = solveApproximationForMultigridUpdated(p,e,t,numOfTriangles(n),k,edge,n,storingA,storingEdge,storeNodeNums,storeHeights);
    storingA{n}                = globalA;
    storeNodeNums{n}           = numOfNodes;
    storeHeights{n}            = height;
    sizeOfDimForMG(n)          = height;
    errorConvergenceRatesMG(n) = MGErrorConvergenceRate;
    numOfMGIterations(n)       = numOfMGIterations1;
    
end

% With Multigrid (Test MG) . . .
MeshLevel                 = zeros(iterations,1);
SizeOfDimForMG            = zeros(iterations,1);
NumOfMGIterations         = zeros(iterations,1);
ErrorConvergenceRateForMG = zeros(iterations,1);

for j = 1:iterations
    MeshLevel(j)                 = meshLevel(j);
    SizeOfDimForMG(j)            = sizeOfDimForMG(j);
    NumOfMGIterations(j)         = numOfMGIterations(j);
    ErrorConvergenceRateForMG(j) = errorConvergenceRatesMG(j);
end

save('MG_Results','MeshLevel','SizeOfDimForMG','NumOfMGIterations','ErrorConvergenceRateForMG');

table(MeshLevel,SizeOfDimForMG,NumOfMGIterations,ErrorConvergenceRateForMG)