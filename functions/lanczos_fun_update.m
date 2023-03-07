function [Xm, iter, lucky] = lanczos_fun_update(A, U, B, f, tol, it, debug)
% Return the core factor Xm of the approximation of the update f(A + U * B * U') - f(A) = Um * Xm * Um' computed with Lanczos (avoids store full basis)
%
%---------------INPUT----------------------------------------------------------------------------------------------------------------
%
% A               matrix argument  TODO modify the code to allow handle function in place of A 
% U               low-rank factor of the update
% f               handle function to be applied 
% tol             (optional) threshold for the (heuristic) stopping criterion, default tol = 1e-12
% it              (optional) max number of iterations, default it = 100
% debug           (optional) allows to print either the heuristic or the true error (the second one is expensive), default debug = 0  
%
%---------------OUTPUT--------------------------------------------------------------------------------------------------------------
%
% Xm          	  core factor in the low-rank factorization of the outcome
% iter		  number of iterations needed to converge
% lucky		  flag indicating that a lucky breakdown has occurred
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
d = 2;       
Xstop = {};

for j = 1:it

	% Augment the Krylov spaces
	if j == 1 %~exist('Um', 'var')
		[Um, HA, param_A, lucky] = lanczos_krylov(A,  U);
		Cm = Um(:, 1:end - rk)' * U; %TODO maybe  this computation can be optimized (minor)
		Cm = Cm * B * Cm';  
	else
		[Um, HA, param_A, lucky] = lanczos_krylov(Um, HA, param_A);
	end

	% Retrieve the projection of A into the  Krylov subspaces
	Gm = HA(1:end - rk, :);

	nn = size(Gm, 1);
	if nn > size(Cm, 1)	
		Cm(nn, nn) = 0;
	end
	if herm 
		tGm = Gm + (Cm + Cm')/2; % enhance symmetry
	else
		tGm = Gm + Cm;
	end
	% Copute the core factor of the outcome 	
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
	warning('LANCZOS_FUN_UPDATE:: Reached maximum number of iterations')
end





