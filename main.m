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
for iter = 1:iterations
    [p, e, t] = refinemesh(g, p, e, t, 'regular');  
    
    c           = approximation(p, e, t, f);
    error(iter) = L2Error(p, e, t, u, c);
    if(iter > 1)
        logError(iter-1) = log2(error(iter-1) / error(iter));
    end 
    % pdeplot(p, e, t, 'XYData', c, 'ZData', c)
end

% plot(logError');