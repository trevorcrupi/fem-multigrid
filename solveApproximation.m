function c = solveApproximation(p,e,t,numOfTriangles)

    countForA      = 1;
    countForB      = 1;
    ourF           = @(x,y) (cos(pi.*y).*(-pi^2-1)/pi);
    b_i            = @(basis, x, y) ourF(x,y).*(basis(1).*x + basis(2).*y + basis(3));
    integrandFunct = @(basisI,basisJ,x,y) (basisI(1).*basisJ(1) + basisI(2).*basisJ(2)) + (basisI(1).*x + basisI(2).*y + basisI(3)).*(basisJ(1).*x + basisJ(2).*y + basisJ(3));

    
%   numOfTriangles = size(t,2);
    for i = 1:numOfTriangles
        columnVector = t(1:3, i); % Vector with points of the triangle
        localPhiCoeffs = getLocalPhiCoeffs(p,columnVector);

        % Use triquad to get localA matrix. (Use 8 because Minah said to.)
        % X,Y,Wx,Wy ??????????
        [X,Y,Wx,Wy]    = triquad(8, [p(1,columnVector(1)) p(2,columnVector(1)); p(1,columnVector(2)) p(2,columnVector(2)); p(1,columnVector(3)) p(2,columnVector(3))]); %triquad(8,[solvePhi(1,1) solvePhi(1,2); solvePhi(2,1) solvePhi(2,2); solvePhi(3,1) solvePhi(3,2)]);


        localA = zeros(3,3);
        for k = 1:3
            for j = 1:3
                localA(k,j) = Wx' * integrandFunct(localPhiCoeffs(:,k),localPhiCoeffs(:,j),X,Y) * Wy; %Why X and Y ???????
            end
        end
        localA;

        for m = 1:3
            for n = 1:3                           % Get the I, J, S for sparse-martix "globalA"
                                                  % Think i1j1, i1j2, . . . , i3j3 in localA.
                AI(countForA) = columnVector(m);  % Global node number of V_m - sparse row#
                AJ(countForA) = columnVector(n);  % Global node number of V_n - sparse column#
                AS(countForA) = localA(m,n);      % (m,n) entry of localA
                countForA = countForA +1;
            end
            %for c = 1:3
                BI(countForB) = columnVector(m);     % Global node number of B_m - sparse row#
                BJ(countForB) = 1;                   % Global node number of B_n - sparse column# - ONLY 1 because vector.
                BS(countForB) = Wx' * b_i(localPhiCoeffs(:,m),X,Y) * Wy;
                countForB = countForB +1;
            %end  WHY NOT???
        end
    end

    height  = size(p,2);      % Number-of-points is the height (# of rows) of the matrix and vectors in the equation.
    globalA = sparse(AI,AJ,AS,height,height);
    globalB = sparse(BI,BJ,BS,height,1);
    c = globalA\globalB

end