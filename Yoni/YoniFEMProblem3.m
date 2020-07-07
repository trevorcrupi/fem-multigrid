clear all

iterations     = 8;
meshLevel      = zeros(iterations,1);
numOfTriangles = zeros(iterations,1);
errorVector    = zeros(iterations,1);
rateOfErrors   = zeros(iterations-1,1);

%load PETForYoni/PETForYoni1.mat
load newEle/new_ele1.mat

k              = 1; % Pick any non-zero integer, positive or negative. 1, -1, 2, -1; 10-(-10) For now.

for n = 1:iterations
    
    mystr = ['PETForYoni/PETForYoni' num2str(n) '.mat'];
    load(mystr);
    meshLevel(n)   = n;
    
    numOfTriangles(n) = size(t,2);
    c                 = solveApproximationForProb3(p,e,t,numOfTriangles,k,n);
    errorVector(n)    = getL2ErrorSquaredforProblem3(c,numOfTriangles(n),p,t);
    
    if(n > 1)
        rateOfErrors(n-1) = log2( errorVector(n-1) / errorVector(n) );
    end
    
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
