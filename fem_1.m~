% Load MyDomain file
load MyDomain.mat;
g = decsg(gd, sf, ns);
model = createpde(1);
geometryFromEdges(model, g);
[p, e, t] = initmesh(g, 'hmax', inf);
pdemesh(p, e, t);
 
% Test functions %
u = @(x, y) -(1/pi).*cos(pi.*y);
f = @(x, y) ((-pi^2-1)/pi).*(cos(pi.*y));

triangles    = size(t, 2);
globalAIndex = 1;
globalBindex = 1;
triquadNodes = 8;

% global basis is a slightly confusing variable. Explanation below: %
%{ 
In order to calculate L2 Error, we need to keep track of ALL our local
basis functions, so we can triquad over the sum of C_i*phi_i in our
triangles. I did not want to create an entire for loop just for this, so I 
am creating a column vector of zeros, then I am just concatenating the 
3x3 basis matrix to it. Thus, we ALWAYS have a column of 0's at the
beginning, so if you want to access ALL the basis coefficients in a giant
matrix, you do this:

globalBasis(:, 2:3*triangles+1)

I.e., we start from column #2 (we dont care about the zero column), then 
we select all the way up to the last column, which is 3*(the # of triangles),
since each triangle has 3 basis vectors, then we add 1 bc the first column
is all zeros. So for a mesh with 2 triangles, there will be 7 total
columns: 1 column of zeros, 3 columns for each triangle. The end.
%}
globalBasis  = zeros(3, 1);

for i = 1 : triangles
  col = t(1:3, i);
  
  basisMatrix = [
      p(1, col(1)) p(2, col(1)) 1; 
      p(1, col(2)) p(2, col(2)) 1; 
      p(1, col(3)) p(2, col(3)) 1 
  ];

  e1 = [ 1 0 0 ]';
  e2 = [ 0 1 0 ]';
  e3 = [ 0 0 1 ]';
 
  basis = [ linsolve(basisMatrix, e1) linsolve(basisMatrix, e2) linsolve(basisMatrix, e3) ];
  globalBasis = cat(2, globalBasis, basis);
  
  % Integral value of basis i and basis j  
  [ X, Y, Wx, Wy ] = triquad(triquadNodes, basisMatrix(:, 1:2));
  
  localA = zeros(3, 3);
  
  % Compute local A_ij matrix %
  % basis1 = column vector of basis i coefficients ( basis1(1) = a1, basis1(2) = b1, basis1(3) = c1  ) %
  % basis2 = column vector of basis j coefficients ( basis1(1) = a1, basis1(2) = b1, basis2(3) = c2  ) % 
  A_ij = @(basis1, basis2, x, y) ... 
    (basis1(1).*basis2(1) + basis1(2).*basis2(2)) ... % compute gradient %
    + (basis1(1).*x + basis1(2).*y + basis1(3)) ...   % set up expression for phi_i %
    .*(basis2(1).*x + basis2(2).*y + basis2(3));      % set up expression for phi_j %
  
  % Compute b_i column vector %
  % basis = column vector of basis i coefficients ( basis(1) = a, basis(2) = b, basis(3) = c ) %
  b_i = @(basis, x, y) ...
      f(x,y) ...                                % given function f, defined above triangle for loop%
      .*(basis(1).*x + basis(2).*y + basis(3)); % expression for phi_i: ax_by+c%
  
  % Fill in local A matrix using A_ij formula and triquad %
  for k = 1:3 
      for j = 1:3
          localA(k, j) = Wx' * A_ij(basis(:, k), basis(:, j), X, Y) * Wy;
      end
  end
  
  % build global A matrix from local matrices %
  for m = 1:3
      for n = 1:3
          globalSparseRow(globalAIndex) = col(m);
          globalSparseCol(globalAIndex) = col(n);
          globalSparseEntry(globalAIndex) = localA(m, n);
          globalAIndex = globalAIndex+1;
      end 
      % create b matrix %
      bSparseRow(globalBindex) = col(m);
      bSparseCol(globalBindex) = 1;
      bSparseEntry(globalBindex) = Wx' * b_i(basis(:, m), X, Y) * Wy;
      globalBindex = globalBindex+1;
  end 
end

globalAMatrix = sparse(globalSparseRow, globalSparseCol, globalSparseEntry, size(p, 2), size(p, 2));
globalBMatrix = sparse(bSparseRow, bSparseCol, bSparseEntry, size(p, 2), 1);

% Our coefficients for u_h % 
c = globalAMatrix\globalBMatrix;


