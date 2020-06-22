pdepoly([0 1 1 0], [0 0 1 1]);
save MyDomain gd sf ns;
g = decsg(gd, sf, ns);
model = createpde(1);
geometryFromEdges(model, g);
[p, e, t] = initmesh(g, 'hmax', inf);
%[p, e, t] = refinemesh(g, p, e, t, 'regular');
%[p, e, t] = refinemesh(g, p, e, t, 'regular');
pdemesh(p, e, t);

% test function in B matrix
f = @(x, y) ((pi^2-1)/pi).*(cos(pi.*y));


% Create sparse matrix for global size (number of nodes x number of nodes)
globalAMatrix = sparse(size(p, 2), size(p, 2));

%col of t matrix = number of triangles
numOfTriangles = size(t, 2);

% counter
globalAIndex = 1;
globalBindex = 1;

for i = 1 : numOfTriangles
  % 3 indeces pointing to x values (in p matrix)of vertices of one triangle 
  col = t(1:3, i);
  col
  % col1 = x_val nodes, col2 = y_val, col3 =1's
  nodes = [
      p(1, col(1)) p(2, col(1)) 1; 
      p(1, col(2)) p(2, col(2)) 1; 
      p(1, col(3)) p(2, col(3)) 1; 
  ];

  %canonical bases
  e1 = [1;0;0];
  e2 = [0;1;0];
  e3 = [0;0;1];
    
  % columns are a,b,c constants in local bases corresponding to points 
  % phi = ax +by + c
  % e1's solution yields basis corresponding to first node on nodes matrix
  phiCoeff = [ linsolve(nodes, e1) linsolve(nodes, e2) linsolve(nodes, e3) ];
  
  % matrix with col1 = x_val & col2 = y_val of nodes forming triangle
  % 3x2 matrix to be input for triquad
  quadritureMatrix = [nodes(:, 1:2)];
  
  % Minah said 8?
  % Confused by triquad and integrand
  % Once integrand is integrated, yields components of local matrices
  [ X, Y, Wx, Wy ] = triquad(8, nodes(:, 1:2));
 
  localA = zeros(3, 3);
  for k = 1:3 
      for j = 1:3
          
        % gradient phi_k dot gradient phi_j yields a_k*a_j + b_k*b_j
        graddot = phiCoeff(1,k).*phiCoeff(1,j) + phiCoeff(2,k).*phiCoeff(2,j);
        %(a_k)x+(b_k)y+(c_k)
        phi_k = @(x,y)(phiCoeff(1,k).*x + phiCoeff(2,k).*y +phiCoeff(3,k));
        %(a_j)x+(b_j)y+(c_j)
        phi_j = @(x,y)(phiCoeff(1,j).*x + phiCoeff(2,j).*y +phiCoeff(3,j));
        
        integrand = @(x, y) graddot + phi_k(x,y) .* phi_j(x,y);
        
        localA(k,j)=Wx'*feval(integrand,X,Y)*Wy;
      end
  end
  localA;
  
  % AX = B
  % entries of B matrix integral over triangle( function times phi_i )
  b_i = @(phiCoeff, x, y) f(x,y)*(phiCoeff(1).*x + phiCoeff(2).*y + phiCoeff(3));
  
% the builds sparse matrix of global matrix A
for m = 1:3
      for n = 1:3
          % row number from the three entries in col 
          globalSparseRow(globalAIndex) = col(m);
          % col number from the three entries in col 
          globalSparseCol(globalAIndex) = col(n);
          % go into the local (m,n)th spot for the value at above indeces
          globalSparseEntry(globalAIndex) = localA(m, n);
          globalAIndex = globalAIndex+1;
      end 
      
      % row number from the three entries in col
      bSparseRow(globalBindex) = col(m);
      % col number is always one since column vector
      bSparseCol(globalBindex) = 1;
      % entries of sparse matrix are the b_i results using the a,b,c for the associated m
      bSparseEntry(globalBindex) = Wx' * b_i(phiCoeff(:, m), X, Y) * Wy;
      globalBindex = globalBindex+1;
  end 
end

% globalAMatrix = sparse(globalSparseRow, globalSparseCol, globalSparseEntry, size(p, 2), size(p, 2));
% globalBMatrix = sparse(bSparseRow, bSparseCol, bSparseEntry, size(p, 2), 1);
% 
% u_h = globalAMatrix\globalBMatrix;
% u_h
