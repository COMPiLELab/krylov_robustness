function c = compute_centrality(A,type)
% Compute the centrality measure of the nodes of a given graph/network
%-------------------------INPUT----------------------------------------------
%
% A				adjacency matrix of the graph
% type			string indicating the kind of centrality measure (see below)
%
%----------------------------------------------------------------------------
    if strcmp(type,'exp')
        c = diag(expm(A));
    elseif strcmp(type,'res')
        rho = abs(eigs(A, 1));
		alpha = 1/(2 * rho);
		c = (speye(n) - alpha * A) \ ones(n, 1);
    elseif strcmp(type,'eig') 
        [u,lambda] = eigs(A,1);
        c = abs(u);
    elseif strcmp(type,'deg')
        c = sum(A);
    elseif strcmp(type,'pr')
        alpha = .85;
        D = diag(1./sum(A));
        n = size(A,1);
        Pf = @(x) alpha * A*D* x + (1-alpha)*sum(x)*(1/n)*ones(n,1);
        [u,lambda] = eigs(Pf,n,1);
        c = abs(u);
    else
        fprintf('#### Centrality set to eig ####\n');
        [u,lambda] = eigs(A,1);
        c = abs(u);
    end
end
