function Hes = hessianfcn_exp(X, A, Omega, tol, it)                                      
% Computing the Hessian
	n = size(A, 1);
	Hes = zeros(size(Omega, 1));
	XX = sparse(Omega(:, 1), Omega(:, 2), X(:), n, n);
	XX = XX + XX';
	Atilde = A + XX;
	[Um, Xm, Vm, row, col, ~] = multiple_frechet_eval(Atilde, Omega, @exp, tol, it, inf, false);
	for j = 1:size(Omega, 1)
		h = Omega(j, 1); k = Omega(j, 2);
		for l = j:size(Omega, 1)
			Hes(j, l) = Um{row(h)}(Omega(l, 1), 1:size(Xm{j}, 1)) * Xm{j} * Vm{col(k)}(Omega(l, 2), 1:size(Xm{j}, 2))';
		end
		Hes(j+1:end, j) = Hes(j, j+1:end); 
	end	
	Hes = -2*Hes;
end
