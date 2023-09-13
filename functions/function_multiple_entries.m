function [X, iter] = function_multiple_entries(A, omega, f, tol, it, poles, debug)
% Approximates a subset of the entries of f(A) with a Krylov projection method for every (i,j) in omega
%
%---------------INPUT----------------------------------------------------------------------------------------------------------------------------------------
%
% A               matrix argument  TODO modify the code to allow handle function in place of A and A'
% omega           (k x 2) array containing the indices (i, j) of f(A) to be computed
% f               handle function to be applied 
% tol             (optional) threshold for the (heuristic) stopping criterion, default tol = 1e-12
% it              (optional) max number of iterations, default it = 100
% poles           (optional) poles for generating the rational Krylov subspace, default poles = inf (polynomial Kryl.)
% debug           (optional) allows to print either the heuristic or the true error (the second one is expensive), default debug = 0  
%
%---------------OUTPUT---------------------------------------------------------------------------------------------------------------------------------------
%
% E			  vector containg the sought entries of f(A), i.e., E(h) = f(A)(omega(h, 1), omega(h, 2))
% iter		  number of iterations needed to converge
%
%------------------------------------------------------------------------------------------------------------------------------------------------------------
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

k = size(omega, 1);
notconverged = [1:k];
rk = 1;
lp = length(poles);
AT = A'; % Speed up observed when the transposed matrix is saved
% Define handle functions for retrieving the Krylov subspaces corresponding to a certain index of row or col
I = unique(omega(:, 1), 'stable')'; nI = length(I);
J = unique(omega(:, 2), 'stable')'; nJ = length(J); 
row = @(t) sum( (I == t) .* [1:nI]);
col = @(t) sum( (J == t) .* [1:nJ]);

% Use MATLAB built-in function if it exists, on the compressed matrix
if isequal(f, @exp)
	f = @(M) expm(M);
elseif isequal(f, @sin)
	f = @(M) funm(M, @sin);
elseif isequal(f, @cos)
	f = @(M) funm(M, @cos);
elseif isequal(f, @log)
	f = @(M) logm(M);
elseif isequal(f, @sqrt)
	f = @(M) sqrtm(M);
else
	f = @(M) funm(M, f);
end

% lag parameter and buffer of old iterates for the stopping criterion 
d = 3;  
Xstop = {{}};   
for h = 1:k-1  
	Xstop = {Xstop{:}, {}};
end

Um = {};
Vm = {};
HA = {};
KA = {};
HB = {};
KB = {};
Xm = {};

Uaux = zeros(nI, 1); 
Vaux = zeros(nJ, 1);

if debug == 2
	fA = f(full(A));
end

for j = 1:it
	% Augment the Krylov spaces for rows
	for h = I
		if j == 1
			U = zeros(size(A, 1), 1); U(h) = 1; 
		if isscalar(poles) && poles == inf			
			[Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A{row(h)}] = arnoldi_krylov(A,  U);
		else
			error('FUNCTION_MULTIPLE_ENTRIES::Unsupported rational Krylov yet')
		end
			temp = Um{row(h)}' * U;
			Uaux(row(h)) = temp(1); 
		else
			if isscalar(poles) && poles == inf			
				[Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A{row(h)}] = arnoldi_krylov(Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A{row(h)});
			else
			error('FUNCTION_MULTIPLE_ENTRIES::Unsupported rational Krylov yet')
			end
		end

		% Retrieve the projection of A into the  Krylov subspaces
		if isscalar(poles) && poles == inf			
			Gm{row(h)} = HA{row(h)}(1:end - rk, :); 
		else
			Gm{row(h)} = HA{row(h)}(1:end - rk, :)/KA{row(h)}(1:end - rk, :); 
		end
	end

	stop = 1; % flag for stopping the external loop	
	for h = notconverged
		% Projection of e_i
		Cm = zeros(j * lp, 1); Cm(1, 1) = Uaux(row(omega(h, 1))); 
	
		% Apply f to the projected matrix
		Xm{h} = f(Gm{row(omega(h, 1))});

		% Check the stopping criterion
		if j <= d
			Xstop{h}{j} = Xm{h};
			stop = 0;
		else
			nn = size(Xm{h}, 1);
			if nn > size(Xstop{h}{1})
				Xstop{h}{1}(nn, nn) = 0;
			end
			err = norm((Xm{h} - Xstop{h}{1}) * eye(nn, 1));

			if debug == 1 || debug == 2
				fprintf('It: %d, Edge: (%d, %d), err = %e, size Xm = %d\n', j, omega(h, 1), omega(h, 2), err, size(Xm{col(omega(h, 2))}, 1));
				if false && debug == 2
					approx_ea = Um{row(omega(h, 1))}(omega(h, 2), 1:size(Xm{h}, 1)) * Xm{h}(:, 1) * Uaux(row(omega(h, 1)));
					true_err = abs(fA(omega(h, 1), omega(h, 2)) - approx_ea) /fA(omega(h, 1), omega(h, 2)); 
					fprintf('True err = %e\n', true_err);
					if true_err > 1e-4 && j >10
					keyboard
					end
				end
			end
			if err > tol % if one of the entry of the matrix function has not converged yet then we go on
				stop = 0; 
			else
				notconverged = setdiff(notconverged, h); % remove entry from loop
				I = unique(omega(notconverged, 1), 'stable')'; 
				J = unique(omega(notconverged, 2), 'stable')'; 
			end

			Xstop{h} = [{Xstop{h}{2:d}}, Xm{h}];
		end
	end
	if stop == 1
		break
	end
end 

iter = j;
if iter == it
	warning('FUNCTION_MULTIPLE_ENTRIES:: Reached maximum number of iterations')
end
X = zeros(size(omega, 1), 1);
for j = 1:k
	X(j) = Um{row(omega(j, 1))}(omega(j, 2), 1:size(Xm{j}, 1)) * Xm{j}(:, 1) * Uaux(row(omega(j, 1)));
	if debug == 2
		true_err = abs(fA(omega(j, 1), omega(j, 2)) - X(j)) /fA(omega(h, 1), omega(h, 2)); 
		fprintf('True err = %e\n', true_err);
		if true_err > 1e-4 
			keyboard
		end
	end
end



