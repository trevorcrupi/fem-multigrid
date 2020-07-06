clear all
load PETForYoni/PETForYoni1.mat
load newEle/new_ele1.mat


edge           = getEdgeMatrix(p,t);
countForA      = 1;
countForB      = 1;
integrandFunct = @(basis_i, basis_j, r) [ ( curlRZKfor(basis_i(1),basis_i(2),basis_i(3)) * curlRZKfor(basis_j(1),basis_j(2),basis_j(3)) + basis_i * basis_j) *r]; % *? or .*? for final r...
bIntegrand     = @(basis_j,r) [ ourF() * basis_j * r ];


%ourF = ;
% curlRZKforPhi = @(U_r,U_theta,U_z) [ - ????(U_theta)
%                                      ???(U_r)
%                                      ???(U_theta) + U_theta/r + (k/r)*U_r ];

% curlRZKforPsi = @(U_r,U_theta,U_z) [ -(k/r)*U_z
%                                      ???(U_r) - ???(U_z)
%                                      (k/r)*U_r ];

%k = ;
psiVector = @(basis,r,z) [ (basis(2)/k).*r - (basis(1)/k).*r.*z
                            0
                           (basis(3)/k).*r + (basis(1)/k).*r.^2 ];
                       
phiVector = @(basis,r,z) [ - (basis(3)/k) - (basis(1)/k).*r - (basis(2)/k).*z
                           basis(1).*r + basis(2).*z + basis(3)
                                                               0];
                                                           

                                                           





numOfTriangles = size(t,2);
 for i = 1:numOfTriangles

    columnVector = t(1:3, i); % Vector with points of the triangle
    rowVector = new_ele(i, 1:3); % Vector with points of the triangle
    
    localPhiCoeffs = getLocalPhiCoeffs(p,columnVector);
    localPsiCoeffs = getLocalPsiCoeffs(p,rowVector,edge);
    
    
    
    
    
    
    [X,Y,Wx,Wy]    = triquad(8, [p(1,columnVector(1)) p(2,columnVector(1)); p(1,columnVector(2)) p(2,columnVector(2)); p(1,columnVector(3)) p(2,columnVector(3))]);
    localA = zeros(6,6);
    for k = 1:3
        for j = 1:3
            ourPhi_i = phiVector( localPhiCoeffs(:,k), X, Y );
            ourPhi_j = phiVector( localPhiCoeffs(:,j), X, Y );
            localA(k,j) = Wx' * integrandFunct(ourPhi_i, ourPhi_j, X) * Wy;
        end
    end
    for k = 1:3
        for j = 4:6
            ourPhi_i = phiVector( localPhiCoeffs(:,k)  , X, Y );
            ourPsi_j = psiVector( localPsiCoeffs(:,j-3), X, Y );
            localA(k,j) = Wx' * integrandFunct(ourPhi_i, ourPsi_j, X) * Wy;
        end
    end
    for k = 4:6
        for j = 1:3
            ourPsi_i = psiVector( localPsiCoeffs(:,k-3), X, Y );
            ourPhi_j = phiVector( localPhiCoeffs(:,j)  , X, Y );
            localA(k,j) = Wx' * integrandFunct(ourPsi_i, ourPhi_j, X) * Wy;
        end
    end
    for k = 4:6
        for j = 4:6
            ourPsi_i = psiVector( localPsiCoeffs(:,k-3), X, Y );
            ourPsi_j = psiVector( localPsiCoeffs(:,j-3), X, Y );
            localA(k,j) = Wx' * integrandFunct(ourPsi_i, ourPsi_j, X) * Wy;
        end
    end
    localA;
    
    
    
    
    
    for m = 1:6
        for n = 1:6                           % Get the I, J, S for sparse-martix "globalA"
                                              % Think i1j1, i1j2, . . . , i3j3 in localA.
            AI(countForA) = columnVector(m);  % Global node number of V_m - sparse row#
            AJ(countForA) = columnVector(n);  % Global node number of V_n - sparse column#
            AS(countForA) = localA(m,n);      % (m,n) entry of localA
            countForA = countForA +1;
        end
        
        
        BI(countForB) = columnVector(m);     % Global node number of B_m - sparse row#
        BJ(countForB) = 1;                   % Global node number of B_n - sparse column# - ONLY 1 because vector.
        
        if m = 1:3
            bFunct_i = phiVector( localPhiCoeffs(:,m), X, Y );
            BS(countForB) = Wx' * bIntegrand(bFunct_i, X) * Wy;
            countForB = countForB +1;
        end
        if m = 4:6
            bFunct_i = psiVector( localPsiCoeffs(:,m-3), X, Y );
            BS(countForB) = Wx' * bIntegrand(bFunct_i, X) * Wy;
            countForB = countForB +1;
        end
 
    end
    
    
    height  = size(p,2) + size(edge,1) ;      % Number-of-points and number-of-edges is the height (# of rows) of the square-matrix and vectors in the equation.
    globalA = sparse(AI,AJ,AS,height,height);
    globalB = sparse(BI,BJ,BS,height,1);
    c       = globalA\globalB;
    
    
    
 end