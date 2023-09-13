function [f, gr] = fun_and_grad_krylov_exp(X, A, Omega, eA, tol, it, debug)
% Evaluates the objective function f(Delta) = trace(expm(A)) - tr(expm(A + Delta)) and 
% its gradient by means of a polynomial Krylov subspace method and the fact that the
% trace of the Frechet derivative in the direction ij corresponds to the ij entry of the 
% matrix exponential. Delta represents the correction which affects the edges in Omega
% with the weights contained in X
%
%-----------------------------------------INPUT-----------------------------------------
%
% A 				Hermitian matrix argument representing the weighted adjacency matrix
% X					vector containing the weight for each edge of the correction 
% Omega				set of edges affected by the correction (Delta)
% eA				entries of expm(A) corresponding to Omega (needed for the Frechet derivative)
% tol				tolerance for the Krylov projection methods
% it				maximum number of iteration for the Krylov projection methods
% debug				allows debugging prints about the accuracy of the sought quantities
%
%----------------------------------------------------------------------------------------
%

	if ~ishermitian(A)
		error('FUN_AND_GRAD_KRYLOV:: matrix A is not Hermitian')
	end
	gr = zeros(size(Omega, 1), 1);
	n = length(A);
	nrmA = normest(A, 1e-2);


	% If the correction is zero
	if sum(abs(X)) == 0
		f = 0;
		gr = -2*eA;
% ------------------DEBUG-----------------------------------------------------------------------------------------------------------
		if debug 
			gr2 = zeros(size(X, 1), 1);
			I = Omega(:, 1);
			J = Omega(:, 2);
			for k = 1 : length(I)    
				i = I(k); j = J(k);
				Eij = zeros(n,n); Eij(i,j) = 1;
				BB = full([A, Eij; zeros(size(A)), A]);
				eB = expm(BB);
				Lfij = eB(1 : n, n + 1 : end);
				gr2(k)= -2 * trace(Lfij);
			end
			if norm(gr - gr2)/norm(gr2) > 1e-5
				[gr, gr2]
				keyboard
			end
			fprintf('Accuracy function eval: %.2e, Accuracy Frechet eval: %.2e, Accuracy Delta: %.2e\n', f, norm(gr - gr2)/norm(gr2), f)
		end
% ----------------END DEBUG---------------------------------------------------------------------------------------------------------
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

	if false && rank(full(B)) < size(B, 1) % hanlde the case where B is rank deficient; it is currently diseabled because it mess up things and seems not needed...
		[Qu, Ru] = qr(U, 0);
		[V, D] = eig(Ru * B * Ru', 'vector');
		ind = find(abs(D) > eps * max(abs(D)));
		B = diag(D(ind));
		U = Qu * V(:, ind);
	end

	[eXm, ~, ~, Um] = fun_update(A, U, B, @exp, tol * exp(nrmA), it, false);
	f = -trace(eXm);

	DeA = Um(Omega(:, 1), :) * eXm * Um(Omega(:, 2), :)';
	DeA = diag(DeA);
	
	gr = -2*(eA + DeA);
% ------------------DEBUG-----------------------------------------------------------------------------------------------------------
	if debug
		AA = full(A + U * B * U');
		XX = expm(AA) - expm(full(A));
		err = norm(XX - full(Um * eXm * Um'))/norm(XX);
		f2 = -trace(XX);
		gr2 = zeros(size(X, 1), 1);
		I = Omega(:, 1);
		J = Omega(:, 2);
		for k = 1 : length(I)    
			i = I(k); j = J(k);
			Eij = zeros(n,n); Eij(i,j) = 1;
			BB = full([AA, Eij; zeros(size(A)), AA]);
			eB = expm(BB);
			Lfij = eB(1 : n, n + 1 : end);
			gr2(k)= -2 * trace(Lfij);
		end
		if norm(gr - gr2)/norm(gr2) > 1e-5
			[gr, gr2]
			keyboard
		end	
		fprintf('Accuracy func. eval: %.2e, Accuracy Frechet eval: %.2e, Accuracy Delta: %.2e\n', abs(f - f2)/abs(f2), norm(gr - gr2)/norm(gr2), err)
	end
% ----------------END DEBUG---------------------------------------------------------------------------------------------------------
end
