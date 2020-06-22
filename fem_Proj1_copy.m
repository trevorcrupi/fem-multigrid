pdepoly([0 1 1 0], [0 0 1 1]);
save MyDomain gd sf ns;
g = decsg(gd, sf, ns);
model = createpde(1);
geometryFromEdges(model, g);
[p, e, t] = initmesh(g, 'hmax', inf);
[p, e, t] = refinemesh(g, p, e, t, 'regular');
%[p, e, t] = refinemesh(g, p, e, t, 'regular');
pdemesh(p, e, t);

% Create sparse matrix for global size (number of nodes x number of nodes)
globalAMatrix = sparse(size(p, 2), size(p, 2));

%col of t matrix = number of triangles
numOfTriangles = size(t, 2);
for i = 1 : numOfTriangles
    
  % col1 = x_val nodes, col2 = y_val, col3 =1's
  nodes = [
      p(1, t(1, i)) p(2, t(1, i)) 1; 
      p(1, t(2, i)) p(2, t(2, i)) 1; 
      p(1, t(3, i)) p(2, t(3, i)) 1; 
  ];
end
  %canonical bases
  e1 = [ 1 0 0 ]';
  e2 = [ 0 1 0 ]';
  e3 = [ 0 0 1 ]';
    
  % columns are a,b,c constants in local bases corresponding to points 
  % phi = ax +by + c
  % e1's solution yields basis corresponding to first node on nodes matrix
  basis = [ linsolve(nodes, e1) linsolve(nodes, e2) linsolve(nodes, e3) ];
  
  % matrix with col1 = x_val & col2 = y_val of nodes forming triangle
  % 3x2 matrix to be input for triquad
  quadritureMatrix = [nodes(:, 1:2)];
  
  % Minah said 8?
  % Confused by triquad and integrand
  % Once integrand is integrated, yields components of local matrices
  [ X, Y, Wx, Wy ] = triquad(8, quadritureMatrix);
 
  localA = zeros(3, 3);
  for i = 1:3 
      for j = 1:3
          
        % gradient phi_i dot gradient phi_j yields a_i*a_j + b_i*b_j
        graddot = basis(1,i)*basis(1,j) + basis(2,i)*basis(2,j);
        %(a_i)x+(b_i)y+(c_i)
        phi_i = @(x,y)(basis(1,i)*x + basis(2,i)*y +basis(3,i));
        %(a_j)x+(b_j)y+(c_j)
        phi_j = @(x,y)(basis(1,j)*x + basis(2,j)*y +basis(3,j));
        
        integrand = @(x, y) graddot + phi_i(x,y) * phi_j(x,y);
        
        localA(i,j)=Wx'*feval(integrand,X,Y)*Wy;
      end
  end
  
 
