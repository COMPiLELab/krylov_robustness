function [Xm, iter, lucky] = lanczos_robustness_update(A, U, B, tol, it, debug)
% Return the approximation of trace(expm(A + U * B * U') - expm(A)) computed with the Lanczos method 
%
%---------------INPUT----------------------------------------------------------------------------------------------------------------
%
% A               Hermitian matrix argument  TODO modify the code to allow handle function in place of A 
% U               low-rank factor of the update
% f               handle function to be applied 
% tol             (optional) threshold for the (heuristic) stopping criterion, default tol = 1e-12
% it              (optional) max number of iterations, default it = 100
% debug           (optional) allows to print either the heuristic or the true error (the second one is expensive), default debug = 0  
%
%---------------OUTPUT--------------------------------------------------------------------------------------------------------------
%
% Xm          Approximation of trace(expm(A + U * B * U') - expm(A))
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

if size(U, 1) <= 130
	fA = full(A);
	fAt = fA + U * B * U';
	fAt = (fAt + fAt')/2;
	Xm = trace(expm(fAt) - expm(fA));
	iter = 0;
	lucky = 0;
	return
end

rk = size(U, 2);
herm = ishermitian(B);

% lag parameter and buffer of old iterates for the stopping criterion 
d = 2;       
Xstop = zeros(1, d);
for j = 1:it

	% Augment the Krylov space
	if j == 1 %~exist('Um', 'var')
		[Um, HA, param_A, lucky] = lanczos_krylov(A,  U);
		Cm = Um(:, 1:end - rk)' * U; %TODO maybe  this computation can be optimized (minor)
		Cm = Cm * B * Cm';  
	else
		[Um, HA, param_A, lucky] = lanczos_krylov(Um, HA, param_A);
	end

	% Retrieve the projection of A into the  Krylov subspace
	Gm = HA(1:end - rk, :);
	nn = size(Gm, 1);
	if nn > size(Cm, 1)	
		Cm(nn, nn) = 0;
	end
	tGm = Gm + Cm;
	if herm % enhance symmetry
		Gm = (Gm + Gm')/2;
		tGm = (tGm + tGm')/2; 
	end
	% Compute the trace of the core factor of the update
	%d1 = sort(eig(Gm + Cm));
	d1 = sort(eig(tGm));
	d2 = sort(eig(Gm));
	Xm = sum(exp(d1) .* (1 - exp(d2 - d1)));
	
	if debug == 3	
		if ~exist('true_quantity', 'var')
			true_quantity = trace(expm(full(A) + U * B * U') - expm(full(A)));
		end
	fprintf('It: %d, err = %e\n', j, abs(true_quantity - Xm)/abs(true_quantity));
	end
	% Check the stopping criterion
	if j <= d
		Xstop(j) = Xm;
	else
		err = abs(Xm - Xstop(1));
		if debug == 1 || debug == 2
			fprintf('It: %d, err = %e, size Gm = %d\n', j, err, size(Gm, 1));
			if debug == 2
				pause
			end
		end
		if err < tol 
			break
		end
		Xstop = [Xstop(2:d), Xm];
	end
	if lucky 
		warning('LANCZOS_ROBUSTNESS_UPDATE:: Detected lucky breakdown')
		break
	end
end 

iter = j;
if iter == it
	warning('LANCZOS_ROBUSTNESS_UPDATE:: Reached maximum number of iterations')
end





