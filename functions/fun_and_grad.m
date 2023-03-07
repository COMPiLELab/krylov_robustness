function [f, gr] = fun_and_grad(X, A, Omega)
% Evaluates the objective function f(Delta) = trace(expm(A)) - tr(expm(A + Delta)) and 
% its gradient with dense arithmetic. Delta represents the 
% correction which affects the edges in Omega with the weights contained in X

	n = length(A);
	D = sparse(Omega(:, 1), Omega(:, 2), X(:)); D(n, n) = 0; D = D + D';
	f = trace(expm(A)) - trace(expm(full(A + D)));


	%% Omega represents a sparse matrix whose nonzeros are the set of entries we 
	%% are allowed to modify
	gr = zeros(size(X, 1), 1);

	I = Omega(:, 1);
	J = Omega(:, 2);
	for k = 1 : length(I)    
		i = I(k); j = J(k);
		Eij = zeros(n,n); Eij(i,j) = 1;
		B = [A + D, Eij; zeros(size(A)), A + D];
		eB = expm(B);
		Lfij = eB(1 : n, n + 1 : end);
		gr(k)= -2 * trace(Lfij);
	end

end



