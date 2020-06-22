% Load MyDomain file
load MyDomain.mat;
g = decsg(gd, sf, ns);
model = createpde(1);
geometryFromEdges(model, g);
[p, e, t] = initmesh(g, 'hmax', inf);
pdemesh(p, e, t);
 
f = @(x, y) ((pi^2-1)/pi).*(cos(pi.*y));


triangles    = size(t, 2);
globalAIndex = 1;
globalBindex = 1;

triquadNodes = 8;

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
  
  % Integral value of basis i and basis j  
  [ X, Y, Wx, Wy ] = triquad(triquadNodes, basisMatrix(:, 1:2));
  
  localA = zeros(3, 3);
  A_ij = @(basis1, basis2, x, y) (basis1(1).*basis2(1) + basis1(2).*basis2(2)) + (basis1(1).*x + basis1(2).*y + basis1(3)).*(basis2(1).*x + basis2(2).*y + basis2(3));
  b_i = @(basis, x, y) f(x,y)*(basis(1).*x + basis(2).*y + basis(3));
  
  for k = 1:3 
      for j = 1:3
          localA(k, j) = Wx' * A_ij(basis(:, k), basis(:, j), X, Y) * Wy;
      end 
  end
  
  for m = 1:3
      for n = 1:3
          globalSparseRow(globalAIndex) = col(m);
          globalSparseCol(globalAIndex) = col(n);
          globalSparseEntry(globalAIndex) = localA(m, n);
          globalAIndex = globalAIndex+1;
      end 
  end 
  
  for q=1:3
    bSparseRow(globalBindex) = col(q);
    bSparseCol(globalBindex) = 1;
    bSparseEntry(globalBindex) = Wx' * b_i(basis(:, q), X, Y) * Wy;
    globalBindex = globalBindex+1;
  end  
end

globalAMatrix = sparse(globalSparseRow, globalSparseCol, globalSparseEntry, size(p, 2), size(p, 2));
globalBMatrix = sparse(bSparseRow, bSparseCol, bSparseEntry, size(p, 2), 1);

u_h = globalAMatrix\globalBMatrix;
u_h