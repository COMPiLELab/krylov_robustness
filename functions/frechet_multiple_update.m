function [Um, Xm, Vm, row, col, iter] = frechet_multiple_update(A, omega, f, tol, it, poles, debug)
% Approximates the Frechet derivative Df(A)(e_i * e_j') = Um * Xm * Vm' with a Krylov projection method for every (i,j) in omega
%
%---------------INPUT----------------------------------------------------------------------------------------------------------------------------------------
%
% A               matrix argument  TODO modify the code to allow handle function in place of A and A'
% omega           (k x 2) array containing the indices (i, j) of the Frechet derivative to be computed
% f               handle function to be applied 
% tol             (optional) threshold for the (heuristic) stopping criterion, default tol = 1e-12
% it              (optional) max number of iterations, default it = 100
% poles           (optional) poles for generating the rational Krylov subspace, default poles = inf (polynomial Kryl.)
% debug           (optional) allows to print either the heuristic or the true error (the second one is expensive), default debug = 0  
%
%---------------OUTPUT---------------------------------------------------------------------------------------------------------------------------------------
%
% row, col        handle functions such that Um{row(omega(h, 1))} and Vm{col(omega(h, 2))} contain the row and column bases corresponding to Xm{h}  
% Um, Xm, Vm      cell arrays of dimension nI,k and nJ respectively, containing the low-rank factorizations of the Frechet derivative in the points in omega
%                 in particular we have:  
%		  if 
%			omega(h, :) = (i, j)  
%		  then
%			Df(A)(e_i * e_j') = Um{row(i)}(:, 1:size(Xm{h}, 1)) * Xm{h} * Vm{col(j)}(:, 1:size(Xm{h}, 2))'
%
%
% Note that Um{row(i)} (and Vm{col(j)} resp.) can contain more col than the number needed for Xm{h}
%
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

if ~exist('compression', 'var')
	compression = 0;
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

Uaux = zeros(nI, 1); % For some strange reasons the starting vector is sometimes represented as either e_1 or -e_1 in the Krylov space, so I need to keep track of this information
Vaux = zeros(nJ, 1);

for j = 1:it
	% Augment the Krylov spaces for rows
	for h = I
		if j == 1
			U = zeros(size(A, 1), 1); U(h) = 1; 
		if isscalar(poles) && poles == inf			
			[Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A{row(h)}] = poly_krylov(A,  U);
		else
			[Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A] = rat_krylov(A,  U, poles);
		end
			temp = Um{row(h)}' * U;
			Uaux(row(h)) = temp(1); 
		else
			if isscalar(poles) && poles == inf			
				[Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A{row(h)}] = poly_krylov(Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A{row(h)});
			else
				[Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A] = rat_krylov(A,  Um{row(h)}, KA{row(h)}, HA{row(h)}, poles, param_A);
			end
		end

		% Retrieve the projection of A into the  Krylov subspaces
		if isscalar(poles) && poles == inf			
			Gm{row(h)} = HA{row(h)}(1:end - rk, :); 
		else
			Gm{row(h)} = HA{row(h)}(1:end - rk, :)/KA{row(h)}(1:end - rk, :); 
		end
	end

	% Augment the Krylov spaces for cols
	for h = J
		if j == 1
			V = zeros(size(A, 2), 1); V(h) = 1; 
			[Vm{col(h)}, KB{col(h)}, HB{col(h)}, param_B] = rat_krylov(AT, V, poles);
			temp = Vm{col(h)}' * V;
			Vaux(col(h)) = temp(1); 
		else
        		[Vm{col(h)}, KB{col(h)}, HB{col(h)}, param_B] = rat_krylov(AT, Vm{col(h)}, KB{col(h)}, HB{col(h)}, poles, param_B);
		end
		% Retrieve the projection of A into the  Krylov subspaces
		Hm{col(h)} = HB{col(h)}(1:end - rk, :)/KB{col(h)}(1:end - rk, :);
	end

	stop = 1; % flag for stopping the external loop	
	for h = notconverged
		% Projection of e_i * e_j'
		Cm = zeros(j * lp); Cm(1, 1) = Uaux(row(omega(h, 1))) * Vaux(col(omega(h, 2))); 
	
		% Build the block upper triangular matrix and apply f to it
		Fm = [Gm{row(omega(h, 1))}, Cm; zeros(size(Hm{col(omega(h, 2))}, 2), size(Gm{row(omega(h, 1))}, 2)), Hm{col(omega(h, 2))}'];
		Fm = f(Fm);

		% Extract core factor of the outcome 	
		Xm{h} = Fm(1:size(Gm{row(omega(h, 1))}, 1), size(Gm{row(omega(h, 1))}, 2) + 1:end);

		% Check the stopping criterion
		if j <= d
			Xstop{h}{j} = Xm{h};
			stop = 0;
		else
			nn = size(Xm{h}, 1);
			if nn > size(Xstop{h}{1})
				Xstop{h}{1}(nn, nn) = 0;
			end
			err = norm(Xm{h} - Xstop{h}{1});
			%err = norm(Xm{h} - blkdiag(Xstop{h}{1}, zeros(size(Xm{h}, 1) - size(Xstop{h}{1}, 1))));
			if debug == 1 || debug == 2
				fprintf('It: %d, Edge: (%d, %d), err = %e, size Hm = %d\n', j, omega(h, 1), omega(h, 2), err, size(Hm{col(omega(h, 2))}, 1));
				if debug == 2
					pause
				end
			elseif debug == 3 || debug == 4	
				C = zeros(size(A)); C(omega(h, 1), omega(h, 2)) = 1;
				Df = f([A C; zeros(size(A)) A]);
				Df = Df(1:size(A, 1), size(A, 2) + 1:end);
				fprintf('It: %d, Edge: (%d, %d), True err = %e, size Hm = %d\n', j, omega(h, 1), omega(h, 2),  norm( Df - Um{row(omega(h, 1))}(:, 1:end - rk) * Xm{h} * Vm{col(omega(h, 2))}(:, 1:end-rk)'), size(Hm{col(omega(h, 2))}, 1));
				if debug == 4
					pause
				end
			end
			if err > tol % if one of the entry of the Frechet derivative has not converged yet then we go on
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
	warning('FRECHET_MULTIPLE_UPDATE:: Reached maximum number of iterations')
end
for j = 1:nI
	Um{j} = Um{j}(:, 1:end - rk);
end

for j = 1:nJ
	Vm{j} = Vm{j}(:, 1:end - rk);
end

