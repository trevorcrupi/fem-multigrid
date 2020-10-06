

iterations = 8;

%The new code . . .
[gd sf ns] = meshgeometryLShapeD(1);
g = decsg(gd, sf, ns);
model = createpde(1);
geometryFromEdges(model, g);
[p, e, t] = initmesh(g, 'hmax', inf);

for n = 1:iterations
    
    %The code that was replaced . . .
%     mystr = ['PETForYoni/PETForYoni' num2str(n) '.mat']; %%%%% NOTE: Didn't update to utilize Matlab's PDE Toolbox -- but it's ok.
%     load(mystr);                                         %%%%% NOTE: Didn't update to utilize Matlab's PDE Toolbox -- but it's ok.
    
    
    %The new code . . .
    if n > 1
        [p, e, t] = refinemesh(g, p, e, t, 'regular');
    end
    
    
    
    edge           = getEdgeMatrix(p,t);
    
    numOfNodes = size(p,2);
    numOfEdges = size(edge,1);
    ZMatrixDraft = zeros(numOfNodes,8)-1;
    
    % Firstly, set a counter into the first column of Z.
    for j = 1:numOfNodes
        ZMatrixDraft(j,1) = 3; % Because first column is counter and second column is node-number.
    end
    
    % Secondly, input the node-numbers into the second column of Z.
    for j = 1:numOfNodes
        ZMatrixDraft(j,2) = j;
    end

    % Then, input edge numbers into Z.
    for i = 1:numOfEdges
        rowForNode1 = edge(i,1);
        rowForNode2 = edge(i,2);
        
        row1Counter = ZMatrixDraft(rowForNode1,1);
        ZMatrixDraft(rowForNode1,row1Counter) = i +numOfNodes; % Want column to be next space with a -1 in it.
        ZMatrixDraft(rowForNode1,1) = ZMatrixDraft(rowForNode1,1) +1;
        
        row2Counter = ZMatrixDraft(rowForNode2,1);
        ZMatrixDraft(rowForNode2,row2Counter) = i +numOfNodes;
        ZMatrixDraft(rowForNode2,1) = ZMatrixDraft(rowForNode2,1) +1;
    end
    
    ZMatrixDraft(:,1) = [];
    ZMatrix = zeros(numOfNodes,7);
    ZMatrix = ZMatrixDraft;
    save(['ZMatrixL', num2str(n),'.mat'],'ZMatrix');
    
end

