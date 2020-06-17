pdepoly([0 1 1 0], [0 0 1 1], 'fem 1 square');
save MyDomain gd sf ns;
g = decsg(gd, sf, ns);
model = createpde(1);
geometryFromEdges(model, g);
[p, e, t] = initmesh(g, 'hmax', inf);
[p, e, t] = refinemesh(g, p, e, t, 'regular');
[p, e, t] = refinemesh(g, p, e, t, 'regular');
pdemesh(p, e, t);