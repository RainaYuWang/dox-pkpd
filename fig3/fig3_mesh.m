function mesh = fig3_mesh(times, wash, step)
% Include requested times and the washout in the piecewise quadrature mesh.
cuts = unique([times(:);wash;2]);
cuts = cuts(cuts>=times(1) & cuts<=times(end));
edges = cuts(1);
for j = 2:numel(cuts)
    localStep = step;
    if cuts(j)<=2, localStep=step/8; end
    count = ceil((cuts(j)-cuts(j-1))/localStep);
    a = linspace(cuts(j-1),cuts(j),count+1)';
    edges = [edges; a(2:end)]; %#ok<AGROW>
end
[~, mesh.indices] = ismember(times,edges);
mesh.dt = diff(edges);
nodes = (1+[-0.8611363115940526,-0.3399810435848563, ...
    0.3399810435848563,0.8611363115940526])/2;
mesh.nodes = edges(1:end-1) + mesh.dt.*nodes;
mesh.weights = [0.3478548451374538;0.6521451548625461; ...
    0.6521451548625461;0.3478548451374538]/2;
% Integrals of the interpolating basis from the left endpoint to each node.
for k = 1:4
    others = nodes([1:k-1,k+1:4]);
    coeff = poly(others)/prod(nodes(k)-others);
    mesh.partial(k,:) = polyval(polyint(coeff),nodes);
end
end
