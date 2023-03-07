function alpha = nm_line_search(f, grad_curr,  alpha, fmax)
delta = 0.8;
gamma1 = 1e-3;
gamma2 = 1;
while f(alpha) < fmax - (alpha * gamma1 + alpha^2 * gamma2) * norm(grad_curr)^2
	alpha = alpha * delta;
end





