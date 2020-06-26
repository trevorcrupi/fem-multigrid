%{
    Name: L2Error
    Auth: Trevor
    Date: 06/24/2020
    Desc: FUNCTION THAT COMPUTES L2 ERROR BETWEEN u and u_h
    Vars:
        - p: XY node data 
        - e: edge data
        - t: triangle data
        - u: exact function to be approximated
        - c: constants for u_h approximation
    Rets: the L2 Error
    -------------------------------------------------
%}

function [error] = L2Error(p, e, t, u, c)
    error = 0;
    for i = 1 : size(t, 2)
        col = t(1:3, i);

        [~, XYData, phi1, phi2, phi3] = computeLocalBasis(p, col);
        [ X, Y, Wx, Wy ] = triquad(20, XYData);

        u_h = @(x, y) c(col(1)).*phi1(x, y) + c(col(2)).*phi2(x, y) + c(col(3)).*phi3(x, y);

        L2SquaredError = Wx' * (u(X,Y)-u_h(X,Y)).^2 * Wy;
        error = error + L2SquaredError;
    end

    error = sqrt(error);
end