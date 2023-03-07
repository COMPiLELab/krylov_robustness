function [edges, rob_variation, A_new] = krylov_miobi_adaptive(A, k, Q, centrality, order, tol, it, poles, debug, miobi, rescale)
% Make it or break it--break/make edge combined with a Krylov projection scheme for the update 
% of the trace  of the matrix exponential and an adaptive strategy for the choice of the relevant edges
%
%-----------------INPUT----------------------------------------------------------------------------------
%
% A				symmetric adjacency matrix describing the network
% k     		budget of edges to be removed/added
% Q     		(optional) number of edges in the search space, default dmax = max degree of A
% centrality  	ranking to pick top Q edges
% order			(optional) ordering strategy on the set of edges
% tol			(optional) tolerance for the inner Krylov projection scheme, default tol = 1e-12
% it    		(optional) maximum number of iterations for the inner Krylov projection scheme, default it = 100
% poles 		(optional) poles used by the inner Krylov projection scheme, default poles = inf
% debug			(optional) debug option for the inner Krylov projection scheme, default debug = false
% miobi 		(optional) string that allows to switch between 'break' and 'make', default miobi = 'break'
% rescale   (optional) value used for rescaling the adjacency matrix (e.g. A has been rescaled to have norm 1), default rescale = 1
%
%-----------------OUTPUT---------------------------------------------------------------------------------
%
% edges 		(k x 2) array containing the edges to be removed/added
% rob_variation reduction/increase of robustness obtained
% Anew			updated adjacency matrix
%
%--------------------------------------------------------------------------------------------------------

if ~issymmetric(A)
	error('KRYLOV_MIOBI:: Adjacency matrix should be symmetric');
end
if ~exist('tol', 'var')
	tol = 1e-12;
end
if ~exist('it', 'var')
	it = min(100, size(A, 1));
end
if ~exist('debug', 'var')
	debug = 0;
end
if ~exist('poles', 'var')
	poles = inf;
end
if ~exist('Q', 'var') || Q == 0
	Q = max(sum(A,1));
end
if ~exist('order', 'var')
	order = 'mult';
end
if ~exist('debug', 'var')
	debug = false;
end
if ~exist('miobi', 'var')
	miobi = 'break';
end
if strcmp(miobi, 'break') && nnz(A) < 2*k 
	error('KRYLOV_MIOBI:: edges to be removed are more than edges in the network')
end
if ~exist('rescale', 'var')
	rescale = 1;
end

rob_variation = 0;
edges = [];

if strcmp(miobi,'make')
    for j = 1 : k 
        % fprintf('Q=%d\n', Q)
		if j == 1
        	top_edges = find_top_missing_edges(A, centrality, Q + k, order); %%% TODO: try update centrality at each step
		else
			n = length(top_edges);
			ind = find(prod(top_edges == tmp_edges, 2));		% remove the previously selected edge from the search space
			top_edges = top_edges([1 : ind - 1, ind + 1 : n], :); 
		end
		E = top_edges(1:Q, :); 
        [tmp_edges, tmp_rob_variation, A] = krylov_miobi(A, 1, E, tol, it, poles, debug, miobi, rescale);
        edges = [edges; tmp_edges];
        rob_variation = rob_variation + tmp_rob_variation;
    end
elseif strcmp(miobi,'break')
    for j = 1 : k 
		if j == 1
        	top_edges = find_top_edges(A, centrality, Q + k, order);
		else
			n = length(top_edges);
			ind = find(prod(top_edges == tmp_edges, 2));		% remove the previously selected edge from the search space
			top_edges = top_edges([1 : ind - 1, ind + 1 : n], :); 
		end
		E = top_edges(1:Q, :);
        [tmp_edges, tmp_rob_variation, A] = krylov_miobi(A, 1, E, tol, it, poles, debug, miobi, rescale);
        edges = [edges; tmp_edges];
        rob_variation = rob_variation + tmp_rob_variation;
    end
end

A_new = A;

end



