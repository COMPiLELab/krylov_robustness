% Test the performances of the interior point scheme incorporating the Hessian for the tuning, rewire, and addition problems on power grid graphs with f=exp

addpath ../functions

%-------------------------- SETTINGS--------------------------------------------------------------------
tol_param = 1e-6; 
poles = inf; 
it = 100; 	       		% max iteration of the Krylov method
debug = 0;
modifiable_edges = 30;  % Size of F, the set of modifiable edges
search_space = 100;		% Size of the first reduction of the search space, based on centrality measures 
heur_method = 'min';	% ordering used to rank edges when centrality measures of node are computed
total_weight = 10; 		% budget for the total weight variation
f = @exp;
methods = ["tuning", "rewire", "add"];
ndense = 500; % if the size of the graph is smaller than this threshold then dense arithmetic is used to compute the initial entries of expm(A)

%--------------------------------Selection of the graph-------------------------------------------------
Problems = load('../datasets_paper/voltage_adjacencies_average_2.mat');
Countries = fieldnames(Problems);
Countries = {Countries{[13, 5, 15, 19, 17, 4, 9, 11, 7, 1]}}; % Restriction to the subset of countries considered in the paper

data = zeros(length(methods), length(Countries));
data_time = data;
data_it = data;

for i = 1 : length(Countries)

    A = getfield(Problems, Countries{i});
    A = A / max(A(:)); 
    trexp = sum(exp(eig(full(A)))); % since the size is moderate we can compute the trace of the exponential via the eigenvalues

    %---------------------------Computation of the set F (called E in the code)------------------------------
    n = length(A);
    centrality = compute_centrality(A, 'eig');
    nrm = f(normest(A, 1e-2)); % Estimated norm of f(A)
    tol = tol_param * nrm;
    
    tic
    % first reduction of the search space
    E = find_top_edges(A, centrality, search_space, heur_method); % existing edges with top centrality measures
    
    % second (and final) reduction of the search space based on the magnitude of the component in the gradient
    if n < ndense
    	expA = expm(full(A));
    	temp = zeros(size(E, 1), 1);
    	for j = 1:size(E, 1)
	    	h = E(j, 1); k = E(j, 2);
	    	temp(j) = expA(h, k);
	    end
    else
        temp = function_multiple_entries(A, E, f, tol, it, poles, debug);
    end	
    
    % since we are interested in minimizing minus the objective function, we look at the most positive entries in the gradient 
    % same approach for all methods
    [temp, ind] = sort(temp, 'descend');
    E = E(ind(1:modifiable_edges), :);
    eA = temp(1:modifiable_edges);
    t1 = toc;
    %---------------------------------------------------------------------------------------------------------

     fprintf('----\n')
    %---------------------------------------------------------------------------------------------------------
    % Zero Initial guess
    x0 = zeros(modifiable_edges, 1);
    
    b = total_weight;
    for hkl = 1 : length(methods)
        method = methods(hkl);

        fprintf('\n')
        fprintf('Graph: %s\t\t\t\tMethod: %s\n', Countries{i},method);

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
        elseif method == "tuning"
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
			if n < ndense
				expA = expm(full(A));
				temp = zeros(size(E, 1), 1);
				for j = 1:size(E1, 1)
	            	h = E1(j, 1); k = E1(j, 2);
	    			temp(j) = expA(h, k);
            	end
			else
				temp = function_multiple_entries(A, E1, f, tol, it, poles, debug);
			end	

            [temp, ind] = sort(temp, 'descend');
            E1 = E1(ind(1:modifiable_edges/2), :);
            eA = temp(1:modifiable_edges/2);
    		
            if n < ndense
				temp = zeros(size(E2, 1), 1);
				for j = 1:size(E2, 1)
	            	h = E2(j, 1); k = E2(j, 2);
	           		temp(j) = expA(h, k);
            	end
			else
				temp = function_multiple_entries(A, E2, f, tol, it, poles, debug);
			end	

            [temp, ind] = sort(temp, 'descend');
            E2 = E2(ind(1:modifiable_edges/2), :);
            E = [E1;E2];
            eA = [eA; temp(1:modifiable_edges/2)];

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
            if n < ndense
                expA = expm(full(A));
				temp = zeros(size(E, 1), 1);
				for j = 1:size(E, 1)
	            	h = E(j, 1); k = E(j, 2);
	            	temp(j) = expA(h, k);
            	end
			else
				temp = function_multiple_entries(A, E, f, tol, it, poles, debug);
			end	
      
            % since we are interested in minimizing minus the objective function, we look at the most positive entries in the gradient 
            % same approach for all methods
            [temp, ind] = sort(temp, 'descend');
            E = E(ind(1:modifiable_edges), :);
            eA = temp(1:modifiable_edges);
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

		% fmincon options
		maxiter = 200;     % max iteration of LBFGS
		options = optimoptions('fmincon',                       ...
								'SpecifyObjectiveGradient',true, ...
								'Display','off', ... 
								'HessianFcn', @(x, lambda) hessianfcn_exp(x, A, E, tol, it),  ...
								'ObjectiveLimit', -1e100, ...
								'ConstraintTolerance', 1e-10, ...
								'MaxIterations', maxiter);         
        
        tic
        
        [x, fval, exitflag, output, lambda, grad] =  fmincon(        ...
                    @(xx) fun_and_grad_krylov_exp(xx, A, E, eA, tol, it, false),               ... %@(x) fun_and_grad(x,A,Omega),          ...
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
        data(hkl, i) = -fval/trexp;
        data_time(hkl, i) = t1+t2;
        data_it(hkl, i) = output.iterations;
        
        filepath = sprintf('../Results/results_weighted_exp_hessian_%s.csv', string(date));
	    save(filepath,'results','-mat');
        
        fprintf('\n');
        fprintf('%d\t%s\t\t%s\t\t\t\t%.2f%%\t\t%.2f\t It: %d\n', n,Countries{i},method,(-fval/trexp)*100,t1+t2,output.iterations);

	end
end

dlmwrite('../Results/results_weighted_exp_hessian_score.dat', data, '\t');
dlmwrite('../Results/results_weighted_exp_hessian_time.dat', data_time, '\t');
dlmwrite('../Results/results_weighted_exp_hessian_it.dat', data_it, '\t');


