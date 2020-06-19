% Load MyDomain file
load MyDomain.mat;
g = decsg(gd, sf, ns);
model = createpde(1);
geometryFromEdges(model, g);
[p, e, t] = initmesh(g, 'hmax', inf);
pdemesh(p, e, t);
 
f = @(x, y) ((pi^2-1)/pi).*(cos(pi.*y));
numOfTriangles = size(t, 2);
count = 1;
bCount = 1;
for i = 1 : numOfTriangles
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
  quadritureMatrix = [ p(1, col(1)) p(2, col(1)); p(1, col(2)) p(2, col(2)); p(1, col(3)) p(2, col(3)) ];
  
  [ X, Y, Wx, Wy ] = triquad(8, quadritureMatrix);
  
  localAIntegrand = @(basis1, basis2, x, y) (basis1(1).*basis2(1) + basis1(2).*basis2(2)) + (basis1(1).*x + basis1(2).*y + basis1(3)).*(basis2(1).*x + basis2(2).*y + basis2(3));
  
  localA = zeros(3, 3);
  for k = 1:3 
      for j = 1:3
          localA(k, j) = Wx' * localAIntegrand(basis(:, k), basis(:, j), X, Y) * Wy;
      end 
  end
  
  for m = 1:3
      for n = 1:3
          globalNodeRowMatrix(count) = col(m);
          globalNodeColMatrix(count) = col(n);
          entryMatrix(count) = localA(m, n);
          count = count+1;
      end 
  end 
  
  bIntegrand = @(basis, x, y) f(x,y)*(basis(1).*x + basis(2).*y + basis(3));
  for q=1:3
    bRowNode(bCount) = col(q);
    bColumnNode(bCount) = 1;
    bColumnEntry(bCount) = Wx' * bIntegrand(basis(:, q), X, Y) * Wy;
    bCount = bCount+1;
  end  
end

globalAMatrix = sparse(globalNodeRowMatrix, globalNodeColMatrix, entryMatrix, size(p, 2), size(p, 2));
globalBMatrix = sparse(bRowNode, bColumnNode, bColumnEntry, size(p, 2), 1);

u_h = globalAMatrix\globalBMatrix;
