function [tr_new, res, it] = mc_trace(Afun, n, tol, maxit, isAreal, debug)
% Trace estimation of the (possibly black-box) matrix Afun via a Monte Carlo approach (convergence ~ 1/maxit) 
%
%--------------------------------------------INPUT----------------------------------------------
%
% Afun		matrix or handle function for computing matvecs
% n			size of the matrix	
% tol		(optional) relative threshold for detecting convergence, default tol = 1e-3
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
	tol = 1e-3;
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
m = 10;

if debug == 1
		fprintf('------------- Trace estimation convergence history -------------\n');
end
K = ceil(maxit/(3*m));
for it = 1:K
	S = sign(randn(n, m));
	G = sign(randn(n, m));
	[Q, ~] = qr(Afun(S), 0);
	tr = tr + trace(Q' * Afun(Q));
	aux = @(x) x - Q * (Q' * x);
	Afun = @(x) aux(Afun(aux(x)));
	tr_new = tr + trace(G' * Afun(G))/m;
	res = abs(tr_new-tr_old)/max(abs(tr_new), abs(tr_old));
	if debug == 1
		fprintf('Number of quadrature pts: %d, Trace estimate: %1.4e, Error: %e\n', it * 3 * m, tr_new, res);
	end
	if res < tol 
		break
	end
	tr_old = tr_new;
end

if isAreal == 1
	tr_new = real(tr_new);
end

