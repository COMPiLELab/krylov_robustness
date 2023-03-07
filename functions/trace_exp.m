function tr = trace_exp(A)
	r = symrcm(A);
	%r = amd(A);
	A = A(r, r);
	Afun = @(x) expmv(1, A, x, [], 'double');
	%tr = mc_trace_hadam(Afun, size(A, 1), 1e-2, round(log2(10^8/nnz(A))), 1);
	tr = mc_trace(Afun, size(A, 1), 1e-4, 1000, 1);
end
