function [Xm, iter, lucky] = trace_fun_update(A, U, B, tol, it, debug, fun)
% Return the approximation of trace(funm(A + U * B * U') - funm(A)) computed with the Lanczos method 
%
%---------------INPUT----------------------------------------------------------------------------------------------------------------
%
% A               Hermitian matrix argument  TODO modify the code to allow handle function in place of A 
% U               low-rank factor of the update
% f               handle function to be applied 
% tol             (optional) threshold for the (heuristic) stopping criterion, default tol = 1e-12
% it              (optional) max number of iterations, default it = 100
% debug           (optional) allows to print either the heuristic or the true error (the second one is expensive), default debug = 0  
% fun			  scalar function (default exp)
%
%---------------OUTPUT--------------------------------------------------------------------------------------------------------------
%
% Xm          Approximation of trace(funm(A + U * B * U') - funm(A))
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

if ~exist('fun', 'var')
	fun = @exp;
end

if size(U, 1) <= 130
	fA = full(A);
	fAt = fA + U * B * U';
	fAt = (fAt + fAt')/2;
	d1 = sort(eig(fAt));
	d2 = sort(eig(fA));
	if isequal(fun, @exp)
		Xm = sum(exp(d1) .* (1 - exp(d2 - d1)));
	else 
		Xm = sum(fun(d1) - fun(d2));
	end	
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
	d1 = sort(eig(tGm));
	d2 = sort(eig(Gm));
	if isequal(fun, @exp)
		Xm = sum(exp(d1) .* (1 - exp(d2 - d1)));
	else 
		Xm = sum(fun(d1) - fun(d2));
	end	
	
	if debug == 3	
		if ~exist('true_quantity', 'var')
			t1 = sort(eig(full(A) + U * B * U'));
			t2 = sort(eig(full(A)));
			if isequal(fun, @exp)
				true_quantity = sum(exp(t1) .* (1 - exp(t2 - t1)));
			else 
				true_quantity = sum(fun(t1) - fun(t2));
			end	
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
        if debug
		    warning('TRACE_FUN_UPDATE:: Detected lucky breakdown')
        end    
		break
	end
end 

iter = j;
if iter == it
	warning('TRACE_FUN_UPDATE:: Reached maximum number of iterations')
end





