function [tr, res, it] = mc_trace_hadam(Afun, n, tol, maxit, isAreal, debug)
% Trace estimation of the (possibly black-box) matrix Afun via Monte Carlo with Hadamard vectors
%
%--------------------------------------------INPUT----------------------------------------------
%
% Afun		matrix or handle function for computing matvecs
% n		size of the matrix	
% tol		(optional) threshold for detecting convergence, default tol = 1e-8
% maxit		(optional) maximum number of iterations, default maxit = 10
% isAreal	(optional) indicates if Afun (and so its trace) is real, default isAreal = 0
% debug		(optional) plots convergence behaviour, default debug = 0
%
%--------------------------------------------OUTPUT----------------------------------------------
%
% tr 		trace of Afun
% res		residual computed for stopping the loop
% it		number of iterations required
%
%------------------------------------------------------------------------------------------------
if ~exist('tol', 'var')
	tol = 1e-8;
end 
if ~exist('maxit', 'var')
	maxit = 10;
end
if ~exist('isAreal', 'var')
	isAreal = 0;
end 
if ~exist('debug', 'var')
	debug = 0;
end
if isfloat(Afun)
	Afun = @(x) Afun * x;
end
tr = 0; tr_old = 0;
p = 1;
%k = 2^ceil(log2(n));

if debug == 1
		fprintf('------------- Trace estimation convergence history -------------\n');
	end
for it = 1:maxit
	tr = tr * p;
	p = 2^it;
	t = floor(n / p);
	r = mod(n, p);
	if it == 1  % Effcient generation of Hadamard vectors
		H = hadamard(2);
		V = [kron(ones(t, 1), H); H(1:r, :)];
Vcomp=V;
	else
		H = hadamard(p/2);
		V = kron(ones(t, 1), [H; -H]);
		if r <= p/2
			V = [V; H(1:r, :)];
		else
			V = [V; H; -H(1:r-p/2, :)];
		end
		Vcomp = [Vcomp, V];
	end
	if it == 1
		tr = trace((V' * Afun(V))); % TODO check if we gain performing scalar products individually
		tr = tr / p;
		if debug == 1
			fprintf('Number of quadrature pts: %d, Trace estimate: %1.4e\n', p, tr);
		end
	else
		tr = tr +  trace((V' * Afun(V))); 
		tr = tr / p;
		res = abs(tr - tr_old)/abs(tr_old);
		if debug == 1
			fprintf('Number of quadrature pts: %d, Trace estimate: %1.4e, Error: %e\n', p, tr, res);
		end
		if res < tol || it == maxit
			break
		end
	end
	tr_old = tr;
end

if isAreal == 1
	tr = real(tr);
end

