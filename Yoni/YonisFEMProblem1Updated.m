clear all

iterations     = 8;
meshLevel      = zeros(iterations,1);
numOfTriangles = zeros(iterations,1);
errorVector    = zeros(iterations,1);
rateOfErrors   = zeros(iterations-1,1);


%The new code . . .
load MyDomain.mat; % Load MyDomain file
g = decsg(gd, sf, ns);
model = createpde(1);
geometryFromEdges(model, g);
[p, e, t] = initmesh(g, 'hmax', inf);

for n = 1:iterations
    
    %The code that was replaced . . .
%     mystr = ['PETForYoni/PETForYoni' num2str(n) '.mat'];
%     load(mystr);


    %The new code . . .
    if n > 1
        [p, e, t] = refinemesh(g, p, e, t, 'regular');
    end
    
    meshLevel(n)      = n;
    numOfTriangles(n) = size(t,2);
    c                 = solveApproximation(p,e,t,numOfTriangles(n));
    errorVector(n)    = getL2Error(c,numOfTriangles(n),p,t);
    
    if(n > 1)
        rateOfErrors(n-1) = log2(errorVector(n-1) / errorVector(n));
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

