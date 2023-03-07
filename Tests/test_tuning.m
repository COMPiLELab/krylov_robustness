%close all 
%clear all

addpath ../functions
addpath ../functions/rktoolbox/

%-------------------------- SETTINGS--------------------------------------------------------------------
tol = 1e-8; 
poles = inf; 
it = 100; 	       		% max iteration of the Krylov method
debug = 0;
modifiable_edges = 30;  % Size of F, the set of modifiable edges
search_space = 100;		% Size of the first reduction of the search space, based on centrality measures 
heur_method = 'min';	% ordering used to rank edges when centrality measures of node are computed
total_weight = 10; 		% budget for the total weight variation
% up_bound = 2;			% upper bound for the weight of a single edge
f = @exp;

% fmincon options
maxiter = 200;     % max iteration of LBFGS
options = optimoptions('fmincon',                       ...
                       'SpecifyObjectiveGradient',true, ...
                       'Display','off', ... 
                       'HessianApproximation','lbfgs',  ...
                       'ObjectiveLimit', -1e100, ...
                       'ConstraintTolerance', 1e-10, ...
                       'MaxIterations', maxiter); 

% function results = run_tuning_on_datasets(modifiable_edges, search_space, heur_method, total_weight,up_bound)

%--------------------------------Selection of the graph-------------------------------------------------
Problems = load('../voltage_adjacencies_average_2.mat');
Countries = fieldnames(Problems);



for i = 1 : length(Countries)
    try
    A = getfield(Problems, Countries{i});
    A = A / max(A(:)); 
    n = length(A);
    m = nnz(A)/2;
    fprintf('%s \t\t\tn=%d\tm=%d\n', Countries{i},n,m);

    catch
     fprintf('Error on %s',Countries{i})
    end
end



keyboard


for i = 1 : length(Countries)
    try
    A = getfield(Problems, Countries{i});
    A = A / max(A(:)); 
    trexp = trace_exp(A);

    %---------------------------Computation of the set F (called E in the code)------------------------------
    n = length(A);
    centrality = compute_centrality(A, 'eig');
    nrm = f(normest(A)); % Estimated norm of f(A)
    tol = tol * nrm;
    
    tic
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
    t1 = toc;
    %---------------------------------------------------------------------------------------------------------



     fprintf('----\n')
    %---------------------------------------------------------------------------------------------------------
    % Zero Initial guess
    x0 = zeros(modifiable_edges, 1);
    %x0 = rand(modifiable_edges, 1);
    methods = ["tuning", "tuning2", "rewire", "add"];
    b = total_weight;
    for hkl = 1 : length(methods)
        method = methods(hkl);

        fprintf('\n')
        fprintf('Graph: %s\t\t\t\tMethod: %s', Countries{i},method);

        % fmincon SYNTAX
        % fmincon(f, x0, M,b,[],[],LB,UB,[],options):
        % minimizes f subject to the constrains
        %      M*x <= b  (in our case sum(sum(X)) <= b)
        %      LB <= x <= UB
        if method == "tuning"
            LB = zeros(size(x0)); 
            UB = LB;
            for j = 1:size(LB, 1)
	            LB(j) = -0.5*A(E(j, 1), E(j, 2));
	            UB(j) = -2*LB(j) ;
            end
        elseif method == "tuning2"
            up_bound = 2;			
            LB = zeros(size(x0)); 
            UB = LB;
            for j = 1:size(LB, 1)
	            LB(j) = -A(E(j, 1), E(j, 2));
	            UB(j) = -LB(j) ;
            end
        elseif method == "rewire"
            tic
            E1 = find_top_edges(A, centrality, search_space/2, heur_method); % existing edges with top centrality measures
            E2 = find_top_missing_edges(A, centrality, search_space/2, heur_method); % existing edges with top centrality measures
            		
            [Um, Xm, Vm, row, col, ~] = frechet_multiple_update(A, E1, f, tol, it, poles, debug);
            temp = zeros(size(E1, 1), 1);
            for j = 1:size(E1, 1)
	            h = E1(j, 1); k = E1(j, 2);
	            temp(j) = trace( Vm{col(k)}(:, 1:size(Xm{j}, 2))' * Um{row(h)}(:, 1:size(Xm{j}, 1)) * Xm{j} );
            end
            [temp, ind] = sort(temp, 'descend');
            E1 = E1(ind(1:modifiable_edges/2), :);
    		
            [Um, Xm, Vm, row, col, ~] = frechet_multiple_update(A, E2, f, tol, it, poles, debug);
            temp = zeros(size(E2, 1), 1);
            for j = 1:size(E2, 1)
	            h = E2(j, 1); k = E2(j, 2);
	            temp(j) = trace( Vm{col(k)}(:, 1:size(Xm{j}, 2))' * Um{row(h)}(:, 1:size(Xm{j}, 1)) * Xm{j} );
            end
            [temp, ind] = sort(temp, 'descend');
            E2 = E2(ind(1:modifiable_edges/2), :);
            E = [E1;E2];
            t1=toc;
            
            LB = zeros(size(x0)); 
            UB = LB;
            for j = 1:size(E1, 1)
	            LB(j) = -A(E(j, 1), E(j, 2));
	            UB(j) = -LB(j) ;
            end
            for j = size(E1, 1)+1:size(E1, 1)+size(E2, 1)
	            LB(j) = 0;
	            UB(j) = 1;
            end 

        elseif method == "add"
            tic;
            % first reduction of the search space
            E = find_top_missing_edges(A, centrality, search_space, heur_method); % existing edges with top centrality measures
            
            % second (and final) reduction of the search space based on the magnitude of the component in the gradient
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
            t1 = toc;
            %---------------------------------------------------------------------------------------------------------

            up_bound = 1;			
            LB = zeros(size(x0)); 
            UB = LB;
            for j = 1:size(LB, 1)
	            LB(j) = 0;
	            UB(j) = up_bound;
            end
        end
        
        tic
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
        t2 = toc; 
        results{i,hkl} = {Countries{i},n, method, E,x,A(sub2ind(size(A),E(:,1),E(:,2))), -fval/trexp, t1+t2};
        filepath = sprintf('../Results/results_tuning_%s.csv', string(date));
	    save(filepath,'results','-mat');
        
        fprintf('\n');
        fprintf('%d\t%s\t\t\t\t%s\t\t\t\t%f.2%%\t\t%f.2\n', n,Countries{i},method,(-fval/trexp)*100,t1+t2);
    
    end

     catch
         fprintf('Error on %s',Countries{i})
     end

end



