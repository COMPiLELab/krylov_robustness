function [Um, Xm, row, iter] = sym_frechet_multiple_update(A, Ind, f, tol, it, poles, debug)
% Approximates the Frechet derivative Df(A)(e_i * e_j') = Um * Xm * Vm' with a Krylov projection method for every (i,j) in omega 
% in the case A = A' and omega = Ind x Ind
%
%---------------INPUT----------------------------------------------------------------------------------------------------------------------------------------
%
% A               matrix argument  TODO modify the code to allow handle function in place of A and A'
% Ind             array of length k such that Ind x Ind is the set of indices (i, j) of the Frechet derivative to be computed
% f               handle function to be applied 
% tol             (optional) threshold for the (heuristic) stopping criterion, default tol = 1e-12
% it              (optional) max number of iterations, default it = 100
% poles           (optional) poles for generating the rational Krylov subspace, default poles = inf (polynomial Kryl.)
% debug           (optional) allows to print either the heuristic or the true error (the second one is expensive), default debug = 0  
%
%---------------OUTPUT---------------------------------------------------------------------------------------------------------------------------------------
%
% row	          handle function such that Um{row(i)} and Um{row(j)} contain the row and column bases corresponding to Xm{row(i), row(j)}  
% Um, Xm          cell arrays of dimension k  and k^2 (only k(k+1)/2 non empty) respectively, containing the low-rank factorizations of the Frechet derivative 
%                 in the points in omega; in particular we have:  
%	
%			if (i, j) is in Ind x Ind and i >= j (if j > i exploit symmetry) then
%	 
%			Df(A)(e_i * e_j') = Um{row(i)}(:, size(Xm{row(i), row(j)}, 1)) * Xm{row(i), row(j)} * Um{row(j)}(:, size(Xm{row(i), row(j)}, 2))'
%
%
% Note that Um{row(i)}  can contain more cols than the number needed to be multiplied with  Xm{row(i), row(j)}
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

k = length(Ind);
notconverged = Ind; % Indices for which we have to continue the Arnoldi process (we start with all indices)
notconverged_pairs = zeros(2, k *(k + 1)/2); % Pairs for which we have to continue the approximation of the Kernel factor (we start with all pairs (i, j) in Ind x Ind such that i >= j)
for h = 1:k
	notconverged_pairs(1, h*(h-1)/2+1: h*(h+1)/2) = Ind(h);
	notconverged_pairs(2, h*(h-1)/2+1: h*(h+1)/2) = Ind(1:h);
end
rk = 1;
lp = length(poles);

% Define a handle function for retrieving the Krylov subspaces corresponding to a certain index of row or col
row = @(t) sum( (Ind == t) .* [1:k]);

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
for h = 1:k^2-1  
	Xstop = {Xstop{:}, {}};
end
Xstop = reshape(Xstop, k, k);

Um = {};
HA = {};
KA = {};
Xm = {};

Uaux = zeros(k, 1); % For some strange reasons the starting vector is sometimes represented as either e_1 or -e_1 in the Krylov space, so I need to keep track of this information

for j = 1:it
	% Augment the Krylov spaces 
	for h = notconverged
		if j == 1
			U = zeros(size(A, 1), 1); U(h) = 1; 
			if isscalar(poles) && poles == inf %&& false
				[Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A{row(h)}] = poly_krylov(A,  U);
			else
				[Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A] = rat_krylov(A,  U, poles);
			end
			temp = Um{row(h)}' * U;
			Uaux(row(h)) = temp(1); 
		else
		if isscalar(poles) && poles == inf %&& false
				[Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A{row(h)}] = poly_krylov(Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A{row(h)});
			else
				[Um{row(h)}, KA{row(h)}, HA{row(h)}, param_A] = rat_krylov(A,  Um{row(h)}, KA{row(h)}, HA{row(h)}, poles, param_A);
			end
		end

		% Retrieve the projection of A into the  Krylov subspaces
		if isscalar(poles) && poles == inf %&& false
			Gm{row(h)} = HA{row(h)}(1:end - rk, :);
		else
			Gm{row(h)} = HA{row(h)}(1:end - rk, :)/KA{row(h)}(1:end - rk, :); 
		end
	end

	stop = 1; % flag for stopping the external loop	
	for h = notconverged_pairs
		% Projection of e_i * e_j'
		Cm = zeros(j * lp); Cm(1, 1) = Uaux(row(h(1))) * Uaux(row(h(2))); 
	
		% Build the block upper triangular matrix and apply f to it
		Fm = [Gm{row(h(1))}, Cm; zeros(size(Gm{row(h(2))}, 2), size(Gm{row(h(1))}, 2)), Gm{row(h(2))}'];
		Fm = f(Fm);

		% Extract core factor of the outcome 	
		Xm{h(1), h(2)} = Fm(1:size(Gm{row(h(1))}, 1), size(Gm{row(h(1))}, 2) + 1:end);

		% Check the stopping criterion
		if j <= d
			Xstop{h(1), h(2)}{j} = Xm{h(1), h(2)};
			stop = 0;
		else
			nn = size(Xm{h(1), h(2)}, 1);
			if nn > size(Xstop{h(1), h(2)}{1}, 1)
				Xstop{h(1), h(2)}{1}(nn, nn) = 0;
			end
			err = norm(Xm{h(1), h(2)} - Xstop{h(1), h(2)}{1});
			%err = norm(Xm{h(1), h(2)} - blkdiag(Xstop{h(1), h(2)}{1}, zeros(size(Xm{h(1), h(2)}, 1) - size(Xstop{h(1), h(2)}{1}, 1))));
			if debug == 1 || debug == 2
				fprintf('It: %d, Edge: (%d, %d), err = %e, size Gm = %d\n', j, h(1), h(2), err, size(Gm{row(h(2))}, 1));
				if debug == 2
					pause
				end
			elseif debug == 3 || debug == 4	
				C = zeros(size(A)); C(h(1), h(2)) = 1;
				Df = f([A C; zeros(size(A)) A]);
				Df = Df(1:size(A, 1), size(A, 2) + 1:end);
				fprintf('It: %d, Edge: (%d, %d), True err = %e, size Gm = %d\n', j, h(1), h(2),  norm( Df - Um{row(h(1))}(:, 1:end - rk) * Xm{h(1), h(2)} * Um{row(h(2))}(:, 1:end-rk)'), size(Gm{row(h(2))}, 1));
				if debug == 4
					pause
				end
			end
			if err > tol % if one of the entry of the Frechet derivative has not converged yet then we go on
				stop = 0; 
			else
				notconverged_pairs = setdiff(notconverged_pairs', h', 'rows')'; 	     % remove pair from the loop
				notconverged = unique([notconverged_pairs(1, :), notconverged_pairs(2, :)]); % if needed remove indices for the Krylov spaces generation 
			end

			Xstop{h(1), h(2)} = [{Xstop{h(1), h(2)}{2:d}}, Xm{h(1), h(2)}];
		end
	end
	if stop == 1
		break
	end
end 
iter = j;
if iter == it
	warning('SYM_FRECHET_MULTIPLE_UPDATE:: Reached maximum number of iterations')
end
for j = 1:k
	Um{j} = Um{j}(:, 1:end - rk);
end
if debug == 5
	for h = 1:k
		notconverged_pairs(1, h*(h-1)/2+1: h*(h+1)/2) = Ind(h);
		notconverged_pairs(2, h*(h-1)/2+1: h*(h+1)/2) = Ind(1:h);
	end
	for h = notconverged_pairs
		C = zeros(size(A)); C(h(1), h(2)) = 1;
		Df = f([A C; zeros(size(A)) A]);
		Df = Df(1:size(A, 1), size(A, 2) + 1:end);
		rk = size(Xm{h(1), h(2)}, 2);
		fprintf('It: %d, Edge: (%d, %d), True err = %e\n', iter, h(1), h(2),  norm( Df - Um{row(h(1))}(:, 1:rk) * Xm{h(1), h(2)} * Um{row(h(2))}(:, 1:rk)'));
	end
end

