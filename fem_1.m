pdepoly([0 1 1 0], [0 0 1 1]);
save MyDomain gd sf ns;
g = decsg(gd, sf, ns);
model = createpde(1);
geometryFromEdges(model, g);
[p, e, t] = initmesh(g, 'hmax', inf);
[p, e, t] = refinemesh(g, p, e, t, 'regular');
[p, e, t] = refinemesh(g, p, e, t, 'regular');
pdemesh(p, e, t);

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
  
  basis1 = linsolve(basisMatrix, e1);
  basis2 = linsolve(basisMatrix, e2);
  basis3 = linsolve(basisMatrix, e3);
  
  xIntegralBounds = integralBounds(p(1, col(1)), p(1, col(2)), p(1, col(3)));
  yIntegralBounds = integralBounds(p(2, col(1)), p(2, col(2)), p(2, col(3)));
  
  % Integral value of basis i and basis j
  B1B1_IV = integrate(basis1, basis1, xIntegralBounds, yIntegralBounds);
  B1B2_IV = integrate(basis1, basis2, xIntegralBounds, yIntegralBounds);
  B1B3_IV = integrate(basis1, basis3, xIntegralBounds, yIntegralBounds);
  B2B1_IV = B1B2_IV;
  B2B2_IV = integrate(basis2, basis2, xIntegralBounds, yIntegralBounds);
  B2B3_IV = integrate(basis2, basis3, xIntegralBounds, yIntegralBounds);
  B3B1_IV = B1B3_IV;
  B3B2_IV = B2B3_IV;
  B3B3_IV = integrate(basis3, basis3, xIntegralBounds, yIntegralBounds);
  
  localA = [
    B1B1_IV B1B2_IV B1B3_IV;
    B2B1_IV B2B2_IV B2B3_IV;
    B3B1_IV B3B2_IV B3B3_IV
  ];
end

% Retrieve the proper bounds using some clever if statements
function [bounds] = integralBounds(x1, x2, x3)
    bounds = [x1; x2];
    if x1 > x2
        bounds = [x2; x1];
    end
    
    if x1 == x2
        if x1 > x3
            bounds = [x3; x1];
        else
            bounds = [x1; x3];
        end
    end
end

% Perform the integral over the basis functions 
function [integralValue] = integrate(basis1, basis2, xBounds, yBounds)
    % create variables from basis function 1
    a1 = basis1(1);
    b1 = basis1(2);
    c1 = basis1(3);
    
    % create variables from basis function 2
    a2 = basis2(1);
    b2 = basis2(2);
    c2 = basis2(3);
    
    % get the bounds from the integral (given by integralBounds function)
    x1 = xBounds(1);
    x2 = xBounds(2);
    y1 = yBounds(1);
    y2 = yBounds(2);
    
    % for ease of reading the integralValue expression
    D = a1*a2+b1*b2;
    x = x2-x1;
    y = y2-y1;
       
    % integral that is computed (computed by hand)
    integralValue = (a1*a2 / 3)*(x2^3-x1^3)*y + (a1*b2 / 4)*(x2^2-x1^2)*(y2^2-y1^2) + (a1*c2 / 2)*(x2^2-x1^2)*y + (b1*a2 / 4)*(x2^2-x1^2)*(y2^2-y1^2) + (b1*b2 / 3)*x*(y2^3-y1^3) + (b1*c2 / 2)*x*(y2^2-y1^2) + (a2*c1 / 2)*(x2^2-x1^2)*y + (b2*c1 / 2)*x*(y2^2-y1^2) + (c1*c2)*x*y + D*x*y;
end
