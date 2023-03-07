function nrm = phi_log(a, A, X, difAX, S, EdgesPattern, r, mu)
%---------------------------------------------------------------------------
%   Evaluate the objective function with log-barrier
%
%   a, X		parameters of the modification
%   A 			adjacency matrix of the original graph
%   difAX               starting (centrality) score
%   S			set of nodes where to evaluate the fitting of scores
%   EdgesPattern 	sparsity pattern of the perturbation
%   r			target (centrality) score
%   mu		 	multiplicative parameter of the log-barrier
%---------------------------------------------------------------------------
    edg = sub2ind(size(A),EdgesPattern(:,1),EdgesPattern(:,2));
    %fAX = expm(full(A + X)); v = diag(fAX);
    nrm = norm(r(S) - a .* difAX(S), 2)^2;
    nrm = nrm - mu .* (log(a)+sum(sum(log(A(edg)+X(edg)))));
end
