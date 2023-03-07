function [Um, Xm, iter] = sym_fun_update(A, U, B, f, tol, it, poles, compression, debug)
% Approximates the update f(A + U * B * U') - f(A) = Um * Xm * Um' with a Krylov projection method
%
%---------------INPUT----------------------------------------------------------------------------------------------------------------
%
% A               matrix argument  TODO modify the code to allow handle function in place of A 
% U               low-rank factor of the update
% f               handle function to be applied 
% tol             (optional) threshold for the (heuristic) stopping criterion, default tol = 1e-12
% it              (optional) max number of iterations, default it = 100
% poles           (optional) poles for generating the rational Krylov subspace, default poles = inf (polynomial Kryl.)
% compression     (optional) additional re-compression at the end, default compression = 0
% debug           (optional) allows to print either the heuristic or the true error (the second one is expensive), default debug = 0  
%
%---------------OUTPUT--------------------------------------------------------------------------------------------------------------
%
% Um, Xm          low-rank factorization of the outcome
% iter		  number of iterations needed to converge
%
%-----------------------------------------------------------------------------------------------------------------------------------
if ~exist('tol', 'var')
	tol = 1e-12;
end 

if ~exist('it', 'var')
	it = min(100, size(A, 1));
end

if ~exist('debug', 'var')
	debug = 0;
end

if ~exist('poles', 'var')
	poles = inf;
end

if ~exist('compression', 'var')
	compression = 0;
end

rk = size(U, 2);
herm = ishermitian(B);
% Use MATLAB built-in function if it exists, on the compressed matrix
if isequal(f, @exp)
	f = @(M) expm(M);
elseif isequal(f, @sin)
	f = @(M) sinm(M);
elseif isequal(f, @cos)
	f = @(M) cosm(M);
elseif isequal(f, @log)
	f = @(M) logm(M);
elseif isequal(f, @sqrt)
	f = @(M) sqrtm(M);
else
	f = @(M) funm(M, f);
end

% lag parameter and buffer of old iterates for the stopping criterion 
d = 3;       
Xstop = {};

for j = 1:it

	% Augment the Krylov spaces
	if j == 1 %~exist('Um', 'var')
		if isscalar(poles) && poles == inf 
			[Um, KA, HA, param_A, lucky] = poly_krylov(A,  U);
		else
			[Um, KA, HA, param_A] = rat_krylov(A,  U, poles);
		end
		Cm = Um(:, 1:end - rk)' * U; %TODO maybe  this computation can be optimized (minor)
		Cm = Cm * B * Cm';  
	else
		if isscalar(poles) && poles == inf 
			[Um, KA, HA, param_A, lucky] = poly_krylov(Um, KA, HA, param_A);
		else
			[Um, KA, HA, param_A] = rat_krylov(A,  Um, KA, HA, poles, param_A);
		end
	end

	% Retrieve the projection of A into the  Krylov subspaces
	if isscalar(poles) && poles == inf
		Gm = (HA(1:end - rk, :) + HA(1:end - rk, :)')/2; % enhance symmetry
	else
		Gm = HA(1:end - rk, :)/KA(1:end - rk, :);
		Gm = (Gm + Gm')/2; % enhance symmetry
	end
        %Cm = blkdiag(Cm, zeros(size(Gm, 1) -  size(Cm, 1)));
	nn = size(Gm, 1);
	if nn > size(Cm, 1)	
		Cm(nn, nn) = 0;
	end
	if herm 
		tGm = Gm + (Cm + Cm')/2; % enhance symmetry
	else
		tGm = Gm + Cm;
	end
	% Compute the core factor of the outcome 	
	Xm = f(Gm + Cm) - f(Gm);

	% Check the stopping criterion
	if j <= d
		Xstop = [Xstop, Xm];
	else
		nn = size(Xm, 1);
		Xstop{1}(nn, nn) = 0; 
		err = norm(Xm - Xstop{1});
		%err = norm(Xm - blkdiag(Xstop{1}, zeros(size(Xm, 1) - size(Xstop{1}, 1))));
		if debug == 1 || debug == 2
			fprintf('It: %d, err = %e, size Gm = %d\n', j, err, size(Gm, 1));
			if debug == 2
				pause
			end
		elseif debug == 3 || debug == 4
			if ~exist('Df', 'var')
				Df = f(A + U * B * U') - f(A);
			end
			fprintf('It: %d, True err = %e, size Gm = %d\n', j,  norm(Df - Um(:, 1:end-rk) * Xm * Um(:, 1:end-rk)'), size(Gm, 1));
			if debug == 4
				pause
			end
		end
		if err < tol
			break
		end
		Xstop = [{Xstop{2:d}}, Xm];
	end
	if lucky 
		warning('LANCZOS_FUN_UPDATE:: Detected lucky breakdown')
		break
	end
end 
iter = j;
if iter == it
	warning('SYM_FUN_UPDATE:: Reached maximum number of iterations')
end
if debug == 5
	Df = f(A + U * B * U') - f(A);
	fprintf('It: %d, True err = %e, size Gm = %d\n', j,  norm(Df - Um(:, 1:end-rk) * Xm * Um(:, 1:end-rk)'), size(Gm, 1));
end
Um = Um(:, 1:end - rk);

if compression == 1
	[V, D] = eig(Xm);
	ind = find(diag(abs(S)) > tol);
	Um = Um * W(:, ind);
	Xm = S(ind, ind);
end




