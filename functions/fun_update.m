function [Xm, iter, lucky, Um] = lanczos_fun_update(A, U, B, fun, tol, it, debug)
% Approximate the update fun(A + U * B * U') - fun(A) = Um * Xm * Um'. If only the core factor Xm is required then it uses Lanczos; otherwise (if also the basis Um is returned) the Arnoldi method
%
%---------------INPUT----------------------------------------------------------------------------------------------------------------
%
% A               matrix argument  TODO modify the code to allow handle function in place of A 
% U               low-rank factor of the update
% fun               handle function to be applied 
% tol             (optional) threshold for the (heuristic) stopping criterion, default tol = 1e-12
% it              (optional) max number of iterations, default it = 100
% debug           (optional) allows to print either the heuristic or the true error (the second one is expensive), default debug = 0  
%
%---------------OUTPUT--------------------------------------------------------------------------------------------------------------
%
% Xm          core factor in the low-rank factorization of the outcome
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
if isequal(fun, @exp)
	f = @(M) expm(M);
elseif isequal(fun, @sin)
	f = @(M) funm(M, @sin);
elseif isequal(fun, @cos)
	f = @(M) funm(M, @cos);
elseif isequal(fun, @log)
	f = @(M) logm(M);
elseif isequal(fun, @sqrt)
	f = @(M) sqrtm(M);
else
	f = @(M) funm(M, fun);
	%f = @(M) fun_diag(M, fun);
end


% lag parameter and buffer of old iterates for the stopping criterion 
d = 2;       
Xstop = {};

for j = 1:it

	% Augment the Krylov spaces
	if nargout <= 3
		if j == 1 %~exist('Um', 'var')
			[Um, HA, param_A, lucky] = lanczos_krylov(A,  U);
			Cm = Um(:, 1:end - rk)' * U; %TODO maybe  this computation can be optimized (minor)
			Cm = Cm * B * Cm';  
		else
			[Um, HA, param_A, lucky] = lanczos_krylov(Um, HA, param_A);
		end
	else	
		if j == 1
			[Um, KA, HA, param_A, lucky] = arnoldi_krylov(A,  U);
			Cm = Um(:, 1:end - rk)' * U; %TODO maybe  this computation can be optimized (minor)
			Cm = Cm * B * Cm';  
		else
			[Um, KA, HA, param_A, lucky] = arnoldi_krylov(Um, KA, HA, param_A);
		end
		if size(Um, 2) >= size(Um, 1)/2 % if the Krylov subspace saturate the space then we compute the quantity with dense arithmetic
			Um = eye(size(U, 1));
			Xm = f(full(A + U * B * U')) - f(full(A));
			iter = j;
			return
		end
	end
	% Retrieve the projection of A into the  Krylov subspaces
	Gm = HA(1:end - rk, :);
	Gm = (Gm + Gm')/2; % enhance symmetry
	
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
	Xm = f(tGm) - f(Gm);

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

Um = Um(:, 1:size(Xm, 1));



