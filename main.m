close all
clear all
% Load MyDomain file
load MyDomain.mat;
g = decsg(gd, sf, ns);
model = createpde(1);
geometryFromEdges(model, g);
[p, e, t] = initmesh(g, 'hmax', inf);

% for newVariableLetter = 1:4
%   [p, e, t] = refinemesh(g, p, e, t, 'regular');  
% end 

iterations = 7;

% Test functions %
u = @(x, y) -(1/pi).*cos(pi.*y);
f = @(x, y) ((-pi^2-1)/pi).*(cos(pi.*y));

error    = zeros(1, iterations);
logError = zeros(1, iterations-1);
N        = zeros(1, iterations);
for iter = 1:iterations
    [p, e, t] = refinemesh(g, p, e, t, 'regular');  
    N(iter)   = size(t, 2);
    size(t, 2)
    c           = approximation(p, e, t, f);
    error(iter) = L2Error(p, e, t, u, c);
    if(iter > 1)
        logError(iter-1) = log2(error(iter-1) / error(iter));
    end 
    % pdeplot(p, e, t, 'XYData', c, 'ZData', c)
end

% Log plot of  error over the triangles % 
loglog(N, error);
% hold on
%loglog([10, 160], [10^-3, (10^-3)/16]);