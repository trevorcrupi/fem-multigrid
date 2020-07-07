clear all
load PETForYoni/PETForYoni1.mat
load newEle/new_ele1.mat


edge           = getEdgeMatrix(p,t);
countForA      = 1;
countForB      = 1;
k              = 1; % Pick any non-zero integer, positive or negative. 1, -1, 2, -1; 10-(-10) For now.
%integrandFunct = @(iFunct1, iFunct2, iFunct3, jFunct1, jFunct2, jFunct3, basisICoeffs, basisJCoeffs, r, z, k, idNum) [ getIntegrand(iFunct1, iFunct2, iFunct3, jFunct1, jFunct2, jFunct3, basisICoeffs, basisJCoeffs, r, z, k, idNum) ] ; % [ ( curlRZKfor(basisI) * curlRZKfor(basisJ) + basis_i * basis_j) *r]; % *? or .*? for final r...
U_r            = @(r,z) z - (1/k).*( (r.^3)/3 - (r.^2)/2 );
U_theta        = @(r,z) -k.*z + (r.^3)/3 - (r.^2)/2;
U_z            = @(r,z) r;


ourF_Line1 = @(r,z) k.*(r-1) + U_r(r,z);
ourF_Line2 = @(r,z) -2.*r+1  + U_theta(r,z);
ourF_Line3 = @(r,z) U_z(r,z);

bIntegrand     = @(jFunct1,jFunct2,jFunct3,r,z) ( ourF_Line1(r,z).*jFunct1 + ourF_Line2(r,z).*jFunct2 + ourF_Line1(r,z).*jFunct3 ) .*r ; %(ourF(r,z) * basis_j) * r ];

psiLine1 = @(basis,r,z) [ (basis(2)/k).*r - (basis(1)/k).*r.*z ];
psiLine2 = @() 0;
psiLine3 = @(basis,r) [ (basis(3)/k).*r + (basis(1)/k).*r.^2 ];
                       
phiLine1 = @(basis,r,z) [ - (basis(3)/k) - (basis(1)/k).*r - (basis(2)/k).*z ];
phiLine2 = @(basis,r,z) [ basis(1).*r + basis(2).*z + basis(3) ];
phiLine3 = @() 0;


numOfTriangles = size(t,2);
 for i = 1:numOfTriangles

    columnVector = t(1:3, i); % Vector with points of the triangle
    rowVector = new_ele(i, 1:3); % Vector with points of the triangle
    
    localPhiCoeffs = getLocalPhiCoeffs(p,columnVector);
    localPsiCoeffs = getLocalPsiCoeffs(p,rowVector,edge);
    
    cAndRVector = [columnVector; (rowVector+size(p, 2))'];
    
    % Triquad - Gaussian quadriture is a way to approximate integral.
    [X,Y,Wx,Wy]    = triquad(8, [p(1,columnVector(1)) p(2,columnVector(1)); p(1,columnVector(2)) p(2,columnVector(2)); p(1,columnVector(3)) p(2,columnVector(3))]);
    localA = zeros(6,6);
    for o = 1:3
        for j = 1:3
            ourPhi_i1 = phiLine1( localPhiCoeffs(:,o), X, Y );
            ourPhi_i2 = phiLine2( localPhiCoeffs(:,o), X, Y );
            ourPhi_i3 = phiLine3();
            
            ourPhi_j1 = phiLine1( localPhiCoeffs(:,j), X, Y );
            ourPhi_j2 = phiLine2( localPhiCoeffs(:,j), X, Y );
            ourPhi_j3 = phiLine3();
            localA(o,j) = Wx' * getIntegrand(ourPhi_i1, ourPhi_i2, ourPhi_i3, ourPhi_j1, ourPhi_j2, ourPhi_j3, localPhiCoeffs(:,o), localPhiCoeffs(:,j), X, Y, k, 1) * Wy;
        end
    end
    for o = 1:3
        for j = 4:6
            ourPhi_i1 = phiLine1( localPhiCoeffs(:,o)  , X, Y );
            ourPhi_i2 = phiLine2( localPhiCoeffs(:,o)  , X, Y );
            ourPhi_i3 = phiLine3();
            
            ourPsi_j1 = psiLine1( localPsiCoeffs(:,j-3), X, Y );
            ourPsi_j2 = psiLine2();
            ourPsi_j3 = psiLine3( localPsiCoeffs(:,j-3), X );
            localA(o,j) = Wx' * getIntegrand(ourPhi_i1, ourPhi_i2, ourPhi_i3, ourPsi_j1, ourPsi_j2, ourPsi_j3, localPhiCoeffs(:,o), localPsiCoeffs(:,j-3), X, Y, k, 2) * Wy;
        end
    end
    for o = 4:6
        for j = 1:3
            ourPsi_i1 = psiLine1( localPsiCoeffs(:,o-3), X, Y );
            ourPsi_i2 = psiLine2();
            ourPsi_i3 = psiLine3( localPsiCoeffs(:,o-3), X );
            
            ourPhi_j1 = phiLine1( localPhiCoeffs(:,j)  , X, Y );
            ourPhi_j2 = phiLine2( localPhiCoeffs(:,j)  , X, Y );
            ourPhi_j3 = phiLine3();
            localA(o,j) = Wx' * getIntegrand(ourPsi_i1, ourPsi_i2, ourPsi_i3, ourPhi_j1, ourPhi_j2, ourPhi_j3, localPsiCoeffs(:,o-3), localPhiCoeffs(:,j), X, Y, k, 3) * Wy;
        end
    end
    for o = 4:6
        for j = 4:6
            ourPsi_i1 = psiLine1( localPsiCoeffs(:,o-3), X, Y );
            ourPsi_i2 = psiLine2();
            ourPsi_i3 = psiLine3( localPsiCoeffs(:,o-3), X );
            
            ourPsi_j1 = psiLine1( localPsiCoeffs(:,j-3), X, Y );
            ourPsi_j2 = psiLine2();
            ourPsi_j3 = psiLine3( localPsiCoeffs(:,j-3), X );
            localA(o,j) = Wx' * getIntegrand(ourPsi_i1, ourPsi_i2, ourPsi_i3, ourPsi_j1, ourPsi_j2, ourPsi_j3, localPsiCoeffs(:,o-3), localPsiCoeffs(:,j-3), X, Y, k, 4) * Wy;
        end
    end
    localA;
    
    
    
    
    
    for m = 1:6
        for n = 1:6                           % Get the I, J, S for sparse-martix "globalA"
                                              % Think i1j1, i1j2, . . . , i3j3 in localA.
            AI(countForA) = cAndRVector(m);  % Global node number of V_m - sparse row#
            AJ(countForA) = cAndRVector(n);  % Global node number of V_n - sparse column#
            AS(countForA) = localA(m,n);      % (m,n) entry of localA
            countForA = countForA +1;
        end
        
        BI(countForB) = cAndRVector(m);     % Global node number of B_m - sparse row#
        BJ(countForB) = 1;                   % Global node number of B_n - sparse column# - ONLY 1 because vector.
        
        if m > 0 && m < 4
            bFunct_i1 = phiLine1( localPhiCoeffs(:,m), X, Y );
            bFunct_i2 = phiLine2( localPhiCoeffs(:,m), X, Y );
            bFunct_i3 = phiLine3();
            BS(countForB) = Wx' * bIntegrand(bFunct_i1, bFunct_i2, bFunct_i3, X, Y) * Wy;
        end
        if m > 3 && m < 7
            bFunct_i1 = psiLine1( localPsiCoeffs(:,m-3), X, Y );
            bFunct_i2 = psiLine2();
            bFunct_i3 = psiLine3( localPsiCoeffs(:,m-3), X );
            BS(countForB) = Wx' * bIntegrand(bFunct_i1, bFunct_i2, bFunct_i3, X, Y) * Wy;
        end
        countForB = countForB +1;
 
    end
    
    
    height  = size(p,2) + size(edge,1) ;      % Number-of-points and number-of-edges is the height (# of rows) of the square-matrix and vectors in the equation.
    globalA = sparse(AI,AJ,AS,height,height);
    globalB = sparse(BI,BJ,BS,height,1);
    c       = globalA\globalB;
  
 end