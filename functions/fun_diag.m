function Y = fun_diag(M, f)
	% Apllies the scalar function f to the (hopefully) diagonalizable matrix H, via diagonalization
	
	if ~ishermitian(M)
		if norm(M-M','fro')/norm(M,'fro') < 1e-14
			M = (M + M') * .5;
			[V, D] = eig(full(M), 'vector');
			Y = V * diag(f(D))* V';
		else
			[V, D] = eig(full(M), 'vector');
			Y = V * diag(f(D))/ V;
		end	
	else
			[V, D] = eig(full(M), 'vector');
			Y = V * diag(f(D))* V';	
	end

end
