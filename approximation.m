
%{
    Name: approximation
    Auth: Trevor
    Date: 06/24/2020
    Desc: GIVEN A FUNCTION f, THIS FUNCTION APPROXIMATES A FUNCTION u
    SOLVING LAPLACES EQUATION. USES THE FINITE ELEMENT METHOD.
    Vars:
        - p: XY node data 
        - localNodesL: column vector of global vertices for a triangle (i.e [1 2 4]')
    Rets: A column vector of constants defining the approximate solution in V_h. 
    -------------------------------------------------
%}
function [c] = approximation(p, e, t, f)
    triangles    = size(t, 2);
    AIndex       = 1;
    BIndex       = 1;

    %{ 
        Compute local A_ij matrix 
        --------------------------------------------------------
        basis1 = column vector of basis i coefficients 
            - ( basis1(1) = a1, basis1(2) = b1, basis1(3) = c1  )
        basis2 = column vector of basis j coefficients 
            - ( basis1(1) = a1, basis1(2) = b1, basis2(3) = c2  ) 
    %}
    A_ij = @(basis1, basis2, x, y) ... 
    (basis1(1).*basis2(1) + basis1(2).*basis2(2)) ... % compute gradient %
    + (basis1(1).*x + basis1(2).*y + basis1(3)) ...   % set up expression for phi_i %
    .*(basis2(1).*x + basis2(2).*y + basis2(3));      % set up expression for phi_j %

    %{ 
        Compute b_i column vector 
        ---------------------------------------------
        basis = column vector of basis i coefficients 
        ( basis(1) = a, basis(2) =b, basis(3) = c ) 
    %}
    b_i = @(basis, x, y) ...
        f(x,y) ...                                        % given function f, defined above triangle for loop%
        .*(basis(1).*x + basis(2).*y + basis(3));         % expression for phi_i: ax_by+c%


    for i = 1 : triangles
      nodes = t(1:3, i);

      [basis, XYData] = computeLocalBasis(p, nodes);

      [ X, Y, Wx, Wy ] = triquad(8, XYData);

      % Fill in local A matrix using A_ij formula and triquad %
      localA = zeros(3, 3);
      for k = 1:3 
          for j = 1:3
              localA(k, j) = Wx' * A_ij(basis(:, k), basis(:, j), X, Y) * Wy;
          end
      end

      % build global A matrix from local matrices %
      for m = 1:3
          for n = 1:3
              ARow(AIndex) = nodes(m);
              ACol(AIndex) = nodes(n);
              AEntry(AIndex) = localA(m, n);
              AIndex = AIndex+1;
          end 
          % create b matrix %
          BRow(BIndex) = nodes(m);
          BEntry(BIndex) = Wx' * b_i(basis(:, m), X, Y) * Wy;
          BIndex = BIndex+1;
      end  
    end

    globalAMatrix = sparse(ARow, ACol, AEntry, size(p, 2), size(p, 2));
    globalBMatrix = sparse(BRow, ones(1, size(BRow, 2)), BEntry, size(p, 2), 1);

    % Our coefficients for u_h % 
    c = globalAMatrix\globalBMatrix;
end 