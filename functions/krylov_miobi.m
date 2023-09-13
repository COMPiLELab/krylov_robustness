function [edges, rob, A_new] = krylov_miobi(A, k, E, tol, it, poles, debug, miobi, rescale)
% Make it or break it--break/make edge combined with a Krylov projection scheme for the update 
% of the trace of the matrix exponential
%
%-----------------INPUT----------------------------------------------------------------------------------
%
% A			symmetric adjacency matrix describing the network
% k     	budget of edges to be removed/added
% E     	(optional) subset of edges that can be removed/added, default E = all edges of A
%       	E has to be a (k x 2) array such that E(j, 1) >= E(j, 2) for all j, for exploiting symmetry
% tol		(optional) tolerance for the inner Krylov projection scheme, default tol = 1e-12
% it    	(optional) maximum number of iterations for the inner Krylov projection scheme, default it = 100
% poles 	(optional) poles used by the inner Krylov projection scheme, default poles = inf
% debug		(optional) debug option for the inner Krylov projection scheme, default debug = false
% miobi 	(optional) string that allows to switch between 'break' and 'make', default miobi = 'break'
% rescale   (optional) value used for rescaling the adjacency matrix (e.g. A has been rescaled to have norm 1), default rescale = 1
%
%-----------------OUTPUT---------------------------------------------------------------------------------
%
% edges 	(k x 2) array containing the edges to be removed/added
% rob   	reduction/increase of robustness obtained
% Anew  	updated adjacency matrix
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
if ~exist('E', 'var') || isempty(E)
	[tmp1, tmp2] = find(A);
	ind = find(tmp1 >= tmp2);
	E = [tmp1(ind), tmp2(ind)];
	clear tmp1 tmp2;
else
	%E = [E; [E(:, 2), E(:, 1)]];
	%ind = find(E(:, 1) >= E(:, 2));
	%E = [E(ind, 1), E(ind, 2)];
end
if ~exist('debug', 'var')
	debug = false;
end
if ~exist('miobi', 'var')
	miobi = 'break';
end
if ~exist('rescale', 'var')
	rescale = 1;
end
if strcmp(miobi, 'break') && nnz(A) < 2*k 
	error('KRYLOV_MIOBI:: edges to be removed are more than edges in the network')
end
nE = size(E, 1);
nA = size(A, 1); 
rob = 0;
edges = zeros(0, 2);
chosen = zeros(1, 2); % temporary variable containing the chosen edge to be removed
for j = 1:min(k, nE)
	if strcmp(miobi, 'break')
		mx = [0 inf];
	else
		mx = [0 -inf];
	end
	for h = 1:nE 	
		if E(h, 1) ~= E(h, 2)		% build the low-rank correction to be applied to the argument
			if strcmp(miobi, 'break')
				B = -[0 1; 1 0] / rescale;
			elseif strcmp(miobi, 'make')
				B = [0 1; 1 0] / rescale;
			else
				error('KRYLOV_MIOBI:: not supported option for miobi')
			end	
			U = zeros(nA, 2);
			U(E(h, 1), 1) = 1;		
			U(E(h, 2), 2) = 1;
		else
			if strcmp(miobi, 'break')
				B = -1;
			elseif strcmp(miobi, 'make')
				B = 1;
			else
				error('KRYLOV_MIOBI:: not supported option for miobi')
			end
			U = zeros(nA, 1);
			U(E(h, 1)) = 1;
		end
		tmp = lanczos_robustness_update(A, U, B, tol, it, debug);
%pause
%------------DEBUG-------------------------
%if true
%	tt = trace(expm(full(A)+ U*B*U'))- eA;
%	err = abs(tt -trace(Xm));
%	if err > 1e-3
%		err
%		keyboard
%	end
%end
%------------END-DEBUG-------------------------
		%tmp = trace(Xm); % evaluate variation in the trace 
		if strcmp(miobi, 'break')
			if tmp < mx(2)			% check (partial) maximality of the variation
				mx(1) = h;
				mx(2) = tmp;
				chosen = E(h, :);
			end
		else
			if tmp > mx(2)			% check (partial) maximality of the variation
				mx(1) = h;
				mx(2) = tmp;
				chosen = E(h, :);
			end
		end		
	end

	E = E([1:mx(1)-1, mx(1)+1:end], :);	% update list of available edges
	nE = nE - 1;
	if strcmp(miobi, 'break')
		A(chosen(1), chosen(2)) = 0; 		% remove edge from the network
		A(chosen(2), chosen(1)) = 0;
	else
		A(chosen(1), chosen(2)) = 1; 		% add edge to the network
		A(chosen(2), chosen(1)) = 1;
	end
	edges = [edges; chosen];     		% update list of removed/added edges
	rob = rob + mx(2); 			% update robustness reduction/increase
end

A_new = A;

end



