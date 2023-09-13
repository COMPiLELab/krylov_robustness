function tr = trace_exp(A)
% Compute an estimate of trace(expm(A)) by means of an Hutchinson estimator approach
	%r = symrcm(A);
	%A = A(r, r);
	Afun = @(x) expmv(1, A, x, [], 'double');
	tr = mc_trace(Afun, size(A, 1), 1e-4, 1000, 1);
end
