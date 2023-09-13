function [f, gr] = fun_and_grad_krylov_sinh(X, A, Omega, fun, dfun, dfA, tol, it, debug, fun_M)
% Evaluates the objective function f(Delta) = trace(f(A)) - tr(f(A + Delta)) and 
% its gradient by means of a polynomial Krylov subspace method and the fact that the
% trace of the Frechet derivative in the direction ij corresponds to the ij entry of the 
% derivative matrix function. Delta represents the correction which affects the edges in Omega
% with the weights contained in X
%
%-----------------------------------------INPUT-----------------------------------------
%
% A 				Hermitian matrix argument representing the weighted adjacency matrix
% X					vector containing the weight for each edge of the correction 
% Omega				set of edges affected by the correction (Delta)
% fun				handle function of f
% dfun				handle function of df
% dfA				entries of df(A) corresponding to Omega (needed for the Frechet derivative)
% tol				tolerance for the Krylov projection methods
% it				maximum number of iteration for the Krylov projection methods
% debug				allows debugging prints about the accuracy of the sought quantities
%
%----------------------------------------------------------------------------------------
%
	if ~ishermitian(A)
		error('FUN_AND_GRAD_KRYLOV_FCONNECTIVITY:: matrix A is not Hermitian')
	end
	gr = zeros(size(Omega, 1), 1);
	n = length(A);
	nrmA = normest(A, 1e-2);


	% If the correction is zero
	if sum(abs(X)) == 0
		f = 0;
		gr = -2*dfA;
		return
	end

	% Retrieve a low-rank factorization of X
	aux = unique(Omega(:));
	k = length(aux);
	iaux = sparse(max(aux), 1); % hash table to retrieve the indices in aux
	iaux(aux) = [1:k]';


	U = zeros(n, k);
	B = zeros(k, k);
	for j = 1:k
		U(aux(j), j) = 1;
	end
	for j = 1:size(Omega, 1)
		i1 = full(iaux(Omega(j, 1)));
		i2 = full(iaux(Omega(j, 2)));
		B(i1, i2) = X(j);
		B(i2, i1) = X(j);
	end

	if false && rank(full(B)) < size(B, 1) % handle the case where B is rank deficient; it is currently diseabled because it mess up things and seems not needed...
		[Qu, Ru] = qr(U, 0);
		[V, D] = eig(Ru * B * Ru', 'vector');
		ind = find(abs(D) > eps * max(abs(D)));
		B = diag(D(ind));
		U = Qu * V(:, ind);
	end

	if isequal(fun, @exp)
		[dfXm, ~, ~, Um] = fun_update(A, U, B, @exp, tol * exp(nrmA), it, false);
		f = -trace(dfXm);
	else
		if nargin(dfun) == 1
			[dfXm, ~, ~, Um] = fun_update(A, U, B, dfun, tol * dfun(nrmA), it, false);
		else
			[dfXm, ~, ~, Um] = fun_update(A, U, B, dfun, tol * dfun(nrmA, 0), it, false);
		end	
		[fXm, ~, ~, ~] = fun_update(A, U, B, fun, tol * fun(nrmA), it, false);	
		f = -trace(fXm);
	end

	DdfA = Um(Omega(:, 1), :) * dfXm * Um(Omega(:, 2), :)';
	DdfA = diag(DdfA);
	
	gr = -2*(dfA + DdfA);
end
