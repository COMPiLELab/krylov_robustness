function [Um, Xm, Vm, iter] = fun_update(A, U, V, f, tol, it, poles, compression, debug)
% Approximates the update f(A + UV') - f(A) = Um * Xm * Vm' with a Krylov projection method
%
%---------------INPUT----------------------------------------------------------------------------------------------------------------
%
% A               matrix argument  TODO modify the code to allow handle function in place of A and A'
% U, V            low-rank factorization of the update
% f               handle function to be applied 
% tol             (optional) threshold for the (heuristic) stopping criterion, default tol = 1e-12
% it              (optional) max number of iterations, default it = 100
% poles           (optional) poles for generating the rational Krylov subspace, default poles = inf (polynomial Kryl.)
% compression     (optional) additional re-compression at the end, default compression = 0
% debug           (optional) allows to print either the heuristic or the true error (the second one is expensive), default debug = 0  
%
%---------------OUTPUT--------------------------------------------------------------------------------------------------------------
%
% Um, Xm, Vm      low-rank factorization of the outcome
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

if size(U, 2) ~= size(V, 2)
	error('FUN_UPDATE:: U and V must have the same number of columns')
else
	rk = size(U, 2);
end
AT = A'; % Speed up observed when the transposed matrix is saved
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
	if j == 1
		if isscalar(poles) && poles == inf
			[Um, KA, HA, param_A] = poly_krylov(A,  U);
			[Vm, KB, HB, param_B] = poly_krylov(AT, V);
		else
			[Um, KA, HA, param_A] = rat_krylov(A,  U, poles);
			[Vm, KB, HB, param_B] = rat_krylov(AT, V, poles);
		end
		Cm = Um(:, 1:end - rk)' * U * (Vm(:, 1:end - rk)' * V)'; %TODO maybe  this computation can be optimized (minor)
	else
		if isscalar(poles) && poles == inf
			[Um, KA, HA, param_A] = poly_krylov(Um, KA, HA, param_A);
        		[Vm, KB, HB, param_B] = poly_krylov(Vm, KB, HB, param_B);
		else
			[Um, KA, HA, param_A] = rat_krylov(A,  Um, KA, HA, poles, param_A);
        		[Vm, KB, HB, param_B] = rat_krylov(AT, Vm, KB, HB, poles, param_B);
		end
	end

	% Retrieve the projection of A into the  Krylov subspaces
	if isscalar(poles) && poles == inf
		Gm = HA(1:end - rk, :); 
		Hm = HB(1:end - rk, :);
	else
		Gm = HA(1:end - rk, :)/KA(1:end - rk, :); 
		Hm = HB(1:end - rk, :)/KB(1:end - rk, :);
	end

	% Build the block upper triangular matrix and apply f to it
	Zm = Vm(:, 1:end - rk)' * U * (Vm(:, 1:end - rk)' * V)'; %TODO maybe  this computation can be optimized (minor)
	nn = size(Gm, 1); mm = size(Hm, 1);
	if nn > size(Cm, 1)
		Cm(nn, mm) = 0;
	end
	Fm = [Gm, Cm; zeros(size(Hm, 2), size(Gm, 2)), Hm' + Zm];
	%Fm = [Gm, blkdiag(Cm, zeros(size(Gm, 1) -  size(Cm, 1), size(Hm, 1) - size(Cm, 2))); zeros(size(Hm, 2), size(Gm, 2)), Hm' + Zm];
	Fm = f(Fm);

	% Extract core factor of the outcome 	
	Xm = Fm(1:size(Gm, 1), size(Gm, 2) + 1:end);

	% Check the stopping criterion
	if j <= d
		Xstop = [Xstop, Xm];
	else
		nn = size(Xm, 1);
		Xstop{1}(nn, nn) = 0;
		err = norm(Xm - Xstop{1});
		%err = norm(Xm - blkdiag(Xstop{1}, zeros(size(Xm, 1) - size(Xstop{1}, 1))));
		if debug == 1 || debug == 2
			fprintf('It: %d, err = %e, size Hm = %d\n', j, err, size(Hm, 1));
			if debug == 2
				pause
			end
		elseif debug == 3 || debug == 4
			if ~exist('Df', 'var')
				Df = f(A + U*V') - f(A);
			end
			fprintf('It: %d, True err = %e, size Hm = %d\n', j,  norm(Df - Um(:, 1:end-rk) * Xm * Vm(:, 1:end-rk)'), size(Hm, 1));
			if debug == 4
				pause
			end
		end
		if err < tol
			break
		end
		Xstop = [{Xstop{2:d}}, Xm];
	end
end 
iter = j;
if iter == it
	warning('FUN_UPDATE:: Reached maximum number of iterations')
end
Um = Um(:, 1:end - rk);
Vm = Vm(:, 1:end - rk);

if compression == 1
	[W, S, Z] = svd(Xm);
	r = sum(diag(S) > tol);
	Um = Um * W(:, 1:r);
	Vm = Vm * Z(:, 1:r);
	Xm = S(1:r, 1:r);
end




