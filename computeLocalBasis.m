%{
    Name: computeLocalBasis
    Auth: Trevor
    Date: 06/24/2020
    Desc: COMPUTES THE LOCAL BASIS FOR A TRIANGLE IN A MESH
    Vars:
        - p: XY node data 
        - localNodesL: column vector of global vertices for a triangle (i.e [1 2 4]')
    Rets:
        - basis: a 3x3 matrix, where each column are the basis coefficients
        - XYData: a 3x2 matrix of coordinates for vertices
        - phi: An individual basis function as a handle ready to use
    -------------------------------------------------
%}

function [basis, XYData, phi1, phi2, phi3] = computeLocalBasis(p, localNodes)
    basisMatrix = [
        p(1, localNodes(1)) p(2, localNodes(1)) 1;
        p(1, localNodes(2)) p(2, localNodes(2)) 1;
        p(1, localNodes(3)) p(2, localNodes(3)) 1
    ];

    e1 = [ 1 0 0 ]';
    e2 = [ 0 1 0 ]';
    e3 = [ 0 0 1 ]';

    basis  = [ linsolve(basisMatrix, e1) linsolve(basisMatrix, e2) linsolve(basisMatrix, e3) ];
    phi1   = @(x, y) (basis(1, 1).*x + basis(2, 1).*y + basis(3,1));
    phi2   = @(x, y) (basis(1, 2).*x + basis(2, 2).*y + basis(3,2));
    phi3   = @(x, y) (basis(1, 3).*x + basis(2, 3).*y + basis(3,3));
    XYData = basisMatrix(:, 1:2);
end