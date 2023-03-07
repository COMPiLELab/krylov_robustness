%close all 
%clear all

addpath ../functions

%-------------------------- SETTINGS--------------------------------------------------------------------
dense = false; % check the result with dense computations
tol = 1e-8; 
poles = inf; 
it = 100; 	       		% max iteration of the Krylov method
debug = 0;
modifiable_edges = 30;  % Size of F, the set of modifiable edges
search_space = 100;		% Size of the first reduction of the search space, based on centrality measures 
heur_method = 'min';	% ordering used to rank edges when centrality measures of node are computed
total_weight = 10; 		% budget for the total weight variation
up_bound = 2;			% upper bound for the weight of a single edge
f = @exp;

% fmincon options
maxiter = 200;     % max iteration of LBFGS
options = optimoptions('fmincon',                       ...
                       'SpecifyObjectiveGradient',true, ...
                       'Display','iter-detailed',       ...
                       'HessianApproximation','lbfgs',  ...
                       'ObjectiveLimit', -1e100, ...
                       'ConstraintTolerance', 1e-10, ...
                       'MaxIterations', maxiter); 
%-------------------------------------------------------------------------------------------------------                       
                       
%--------------------------------Selection of the graph-------------------------------------------------
%%%%%% DATASETS  %%%%%%%%%%%%%%
%load ../Datasets/london.mat % in this case comment 'A = Problem.A'
%load ../Datasets/USAir97.mat; 
%load ../Datasets/Journals.mat; 
%load ../Datasets/netscience.mat; 
%load ../Datasets/karate.mat; 
%load ../'MIOBI Codes'/dt_oregon.mat
load('../voltage_adjacencies_average_2.mat');
%% A0 A1 A2 .... A8 are the 9 Oregon matrices of Chan Akoglu Tong
%A = A0; 
A = Italy; 
%A = USA_South;
A = A / max(A(:)); 
%A = Problem.A; 
%--------------------------------------------------------------------------------------------------------

%---------------------------Computation of the set F (called E in the code)------------------------------
n = length(A);
centrality = compute_centrality(A, 'eig');
nrm = f(normest(A)); % Estimated norm of f(A)
tol = tol * nrm;

% first reduction of the search space
E = find_top_edges(A, centrality, search_space, heur_method); % existing edges with top centrality measures

% second (and finale) reduction of the search space based on the magnitude of the component in the gradient
[Um, Xm, Vm, row, col, ~] = frechet_multiple_update(A, E, f, tol, it, poles, debug);
temp = zeros(size(E, 1), 1);
for j = 1:size(E, 1)
	h = E(j, 1); k = E(j, 2);
	temp(j) = trace( Vm{col(k)}(:, 1:size(Xm{j}, 2))' * Um{row(h)}(:, 1:size(Xm{j}, 1)) * Xm{j} );
end
% since we are interested in minimizing minus the objective function, we look at the most positive entries in the gradient 
% same approach for all methods
[temp, ind] = sort(temp, 'descend');
E = E(ind(1:modifiable_edges), :);
%---------------------------------------------------------------------------------------------------------

%---------------------------------------------------------------------------------------------------------

% Zero Initial guess
x0 = zeros(modifiable_edges, 1);
%x0 = rand(modifiable_edges, 1);

% options = optimoptions(@fmincon,'OutputFcn',@outfun);

%% fmincon SYNTAX
%% fmincon(f, x0, M,b,[],[],LB,UB,[],options):
%% minimizes f subject to the constrains
%%      M*x <= b  (in our case sum(sum(X)) <= b)
%%      LB <= x <= UB
b = total_weight;
LB = zeros(size(x0)); 
UB = LB;
for j = 1:size(LB, 1)
	LB(j) = -A(E(j, 1), E(j, 2));
	UB(j) = up_bound + LB(j) ;
end
if dense
[xdense, fval, exitflag, output, lambda, grad] =  fmincon(        ...
            @(xx) fun_and_grad(xx, A, E),               ... %@(x) fun_and_grad(x,A,Omega),          ...
            x0,                                         ...
            ones(size(x0))',                            ...
            b,                                          ...
            [],[],                                      ...
            LB,UB,                                      ...
            [], options);

end
% keyboard
[x, fval, exitflag, output, lambda, grad] =  fmincon(        ...
            @(xx) fun_and_grad_krylov(xx, A, E, tol, it, false),               ... %@(x) fun_and_grad(x,A,Omega),          ...
            x0,                                         ...
            ones(size(x0))',                            ...
            b,                                          ...
            [],[],                                      ...
            LB,UB,                                      ...
            [], options);
XX = full(sparse(E(:, 1), E(:, 2), x(:)));
XX(n, n) = 0; XX = XX + XX';

trexp = trace_exp(A);
 
fprintf('Percentage of approximated Robustness increase: %f\n', -fval/trexp);
if dense
	fprintf('  Edge\t\t Krylov weight\t\t Dense weight\n')
else
	fprintf('  Edge\t\t Krylov weight\t\t Original weight \t\t Weight\n')
end	


for j = 1:size(E, 1)
	if dense	
		fprintf(' (%d, %d) \t    %1.2f \t\t  %1.2f\n', E(j, 1), E(j, 2), x(j), xdense(j))
	else
		fprintf(' (%d, %d) \t    %1.2f \t\t  %1.2f \t\t\t  %1.2f\n', E(j, 1), E(j, 2), x(j), full(A( E(j, 1), E(j, 2)) ), temp(j))
	end
end


