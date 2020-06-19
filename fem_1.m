pdepoly([0 1 1 0], [0 0 1 1]);
save MyDomain gd sf ns;
g = decsg(gd, sf, ns);
model = createpde(1);
geometryFromEdges(model, g);
[p, e, t] = initmesh(g, 'hmax', inf);
% [p, e, t] = refinemesh(g, p, e, t, 'regular');
% [p, e, t] = refinemesh(g, p, e, t, 'regular');
pdemesh(p, e, t);

% Create sparse matrix for global size (number of nodes x number of nodes)
globalAMatrix = sparse(size(p, 2), size(p, 2));

numOfTriangles = size(t, 2);
for i = 1 : numOfTriangles
  col = t(1:3, i);
  
  basisMatrix = [
      p(1, col(1)) p(2, col(1)) 1; 
      p(1, col(2)) p(2, col(2)) 1; 
      p(1, col(3)) p(2, col(3)) 1 
  ];

  e1 = [ 1 0 0 ].';
  e2 = [ 0 1 0 ].';
  e3 = [ 0 0 1 ].';
    
  basis = [ linsolve(basisMatrix, e1) linsolve(basisMatrix, e2) linsolve(basisMatrix, e3) ];
  
  % Integral value of basis i and basis j
  quadritureMatrix = [ 
      p(1, col(1)) p(2, col(1)); 
      p(1, col(2)) p(2, col(2)); 
      p(1, col(3)) p(2, col(3))     
  ];
  
  [ X, Y, Wx, Wy ] = triquad(20, quadritureMatrix);
  
  integrand = @(basis1, basis2, x, y) (basis1(1).*basis2(1) + basis1(2).*basis2(2)) + (basis1(1).*x + basis1(2).*y + basis1(3)).*(basis2(1).*x + basis2(2).*y + basis2(3));
  
  localA = zeros(3, 3);
  for k = 1:3 
      for j = 1:3
          localA(k, j) = Wx' * feval(integrand, basis(:, k), basis(:, j), X, Y) * Wy;
      end 
  end 
  
  localA
end
