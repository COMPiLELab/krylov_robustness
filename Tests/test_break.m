% Test the performances of the various greeedy schemes for the break problem on either benchamrk or transportation graphs (determined by datasets{1} or dataset{2} in the for loop)

clear all

addpath ../functions
addpath ../'MIOBI Codes'/


datasets{1} = dir('../datasets_paper/Misc/*');
datasets{2} = dir('../datasets_paper/Transport/*');
datasets{3} = load('../voltage_adjacencies_average_2.mat');
column_names = {'method', 'dataset', 'n', 'm', ...
                'searchspace_krylovmiobi', 'centrality_order', 'time', 'tr_variation','budget_size'};
Results_TAB = table([],[],[],[],[],[],[],[],[],'VariableNames',column_names);
 
tol = 1e-6; 
poles = inf; 
it = 100; 
debug = 0;
budget_array = floor(linspace(10,100,10)); % budget       
Q_array = [50 250 1000]; % tunes the size of the search space in krylovmiobi_adaptive
	  

debug = false;

for ii = 1:3 % Modify this line to loop over benchmark, transportation graphs, or power grids
    if ii<3 
    	data = datasets{ii}; 
    	l = 3;
    else
    	data = fieldnames(datasets{ii});
    	l = 1;
	end
	for i =  l:length(data) %start from 3 as first two are '.' and '..'
		if ii<3
			filename = sprintf('../datasets_paper/Misc/%s', data(i).name);
			name = split(data(i).name,'.');
			name = string(name{1});
			load(filename); 
			A = Problem.A;
		else
			A = getfield(datasets{ii}, data{i});
			name = string(data{i});
		end	

		%%% make A symmetric and remove selfloops
		A = spones(A + A'); 
		A = A - spdiags(diag(A), 0, size(A, 1), size(A, 2));
		dmax = max(sum(A));
		
		%%% get max connected component
		ind = max_connected_component(A);
		A = A(ind,ind);
		n = size(A,1);
		
		%%% for stopping criterion we estimate norm of exp A
		nrm = exp(normest(A));
	   	trexp = trace_exp(A);

		f = @exp;
		
		%%% We use eigenvector centrality to choose search space
		centrality = compute_centrality(A,'eig');
		centrality_order = {'mult','min'};

		for k = budget_array
		
		    budget_size = k;
			fprintf('Dataset: %s\t n: %d\t budget: %d\t ||exp(A)|| = %.1e\n', name, n, k, nrm);
		    
	 
		    for co = 1:2 % loop on centrality orders
		        %% KRYLOV MIOBI
		        %method = "krylov_miobi";
		        size_search_space_kmiobi = min(1000, nnz(A)/2);

		        %tic;
		        %E = find_top_edges(A, centrality, size_search_space_kmiobi, centrality_order{co});
		        %[krylovedges, tr_variation, A_new] = krylov_miobi(A, k, E, tol * nrm, it, poles, debug,'break');
		        %time_kmiobi= toc;
		       
		        %Results_TAB = [Results_TAB; ...
		        %    {method, name,n,nnz(A)/2,size_search_space_kmiobi,...
		        %    string(centrality_order{co}),time_kmiobi,tr_variation,budget_size}];
				
				% disp(Results_TAB(Results_TAB.dataset == name,:))

		        %edadd = nnz(spones(spones(A) - spones(A_new)));
				if debug 
					fprintf('Number of new edges in KRYLOV MIOBI = %d\n', edadd/2);
				end

		        for Q = Q_array
		            %% ADAPTIVE KRYLOV MIOBI
		            method  = "miobi_adaptive";
		            Q = min(nnz(A)/2 - k, Q);
		            tic;
		            size_search_space_kmiobi_adaptive = Q + k;
		            [edges, delta_trace, A_new] = krylov_miobi_adaptive(A, k, Q, centrality, centrality_order{co}, tol * nrm, it, poles, debug,'break');
		            time_kmiobi_adaptive = toc;
	 
		            Results_TAB = [Results_TAB; ...
		            {method, name,n,nnz(A)/2,size_search_space_kmiobi_adaptive,...
		            string(centrality_order{co}),time_kmiobi_adaptive,delta_trace/trexp,budget_size}];
		        
		            edadd = nnz(spones(spones(A) - spones(A_new)));
					if debug 
						fprintf('Number of new edges in KRYLOV MIOBI ADAPTIVE = %d\n', edadd/2);
					end
		        end
		        
		        %% Random
		        method = "random";
		        tic;
		        E = find_top_edges(A,centrality,size_search_space_kmiobi,centrality_order{co});
		        p = randperm(size(E,1));
		        RANDOMEdges = E(p(1:k),:);
		    	[U, B] = edge2low_rank(RANDOMEdges, n);
		        delta_trace = lanczos_robustness_update(A, full(U), B, tol * nrm);
		        time_random = toc;
		        Results_TAB = [Results_TAB; ...
		            {method, name,n,nnz(A)/2,size_search_space_kmiobi,...
		            string(centrality_order{co}),time_random,delta_trace/trexp,budget_size}];
				if debug 
					fprintf('Random done\n')
				end
		    end
		    
		    %% MIOBI 25 eigs
		    method = "miobi_25";
		    num_eig_miobi = 25; % num of eigenvectors to use in MIOBI
		    tic;
			[~, ~, Aorg] = MIOBIBreakEdge2(A,k,num_eig_miobi);
			[t1, t2, ~] = find(spones(spones(A) - spones(Aorg)));
			[U, B] = edge2low_rank([t1 t2], n);
		    delta_trace = lanczos_robustness_update(A, full(U), B, tol * nrm);
			time_miobi_25 = toc;
		    
		    [~,ind] = sort(centrality,'descend');
		    search_space_miobi = nnz(A)/2;
		    
		    Results_TAB = [Results_TAB; ...
		        {method, name,n,nnz(A)/2,search_space_miobi,...
		        "--",time_miobi_25,delta_trace/trexp,budget_size}];

			if debug 
				fprintf('miobi_25 done\n')
			end
		    
		    %% HEURISRISTIC Arrigo&Benzi,SISC 2016
		    method = "heuristic_mult";
		    tic; 
		    [~,ind] = sort(centrality, 'descend');
		    ind_top_nodes = ind(1:ceil(n/5));
		    Asmall = A(ind_top_nodes,ind_top_nodes); 
		    if nnz(Asmall)<2*k
		        Asmall = A;
		        ind_top_nodes = 1:n;
		    end
		    
		    EE = find_top_edges(Asmall,centrality(ind_top_nodes),k,'mult');
		    HURISTICEdges = [ind_top_nodes(EE(:,1)) ind_top_nodes(EE(:,2))];
		    
		    [U, B] = edge2low_rank(HURISTICEdges, n);
			delta_trace = lanczos_robustness_update(A, full(U), B, tol * nrm);
			time_heuristic = toc;
			Results_TAB = [Results_TAB; ...
		        {method, name,n,nnz(A)/2,k,"mult",time_heuristic,delta_trace/trexp,budget_size}];
			if debug 
				fprintf('heursitic mult done\n')
			end
		    
		    %% HEURISRISTIC
		    method = "heuristic_min";
		    tic; 
		    [~,ind] = sort(centrality, 'descend');
		    ind_top_nodes = ind(1:ceil(n/5));
		    Asmall = A(ind_top_nodes,ind_top_nodes); 
		    if nnz(Asmall) < 2 * k
		        Asmall = A;
		        ind_top_nodes = [1:n]';
		    end
		    
		    EE = find_top_edges(Asmall,centrality(ind_top_nodes),k,'min');
		    HURISTICEdges = [ind_top_nodes(EE(:,1)) ind_top_nodes(EE(:,2))];
		    
		    [U, B] = edge2low_rank(HURISTICEdges, n);
			delta_trace = lanczos_robustness_update(A, full(U), B, tol * nrm);
			time_heuristic = toc;
			Results_TAB = [Results_TAB; ...
		        {method, name,n,nnz(A)/2,k,"min",time_heuristic,delta_trace/trexp,budget_size}];
			if debug 
				fprintf('heuristic min done\n')
			end		   
		    filepath = sprintf('../Results/results_break_%s.csv', string(date));
		    writetable(Results_TAB,filepath)
			disp(Results_TAB(Results_TAB.dataset == name,:))
		    
		end

	end
end

%--------------- Auxiliary functions --------------------
function ind = max_connected_component(A)
    bins = conncomp(graph(A));
    labels = unique(bins);
    maxval = 0; maxindex = 1;
    for i = 1 : length(labels)
        val = sum(bins==labels(i));
        if val > maxval, maxval = val; maxindex = i; end
    end
    ind = (bins==maxindex);
end
%-----------------------------------------------------------
function [U, B] = edge2low_rank(E, n)
% Build the low-rank factorization that removes the edges in E
 	t1 = E(:,1);
    t2 = E(:,2);
    ut1  = unique([t1; t2]);
    U = sparse(ut1, 1:length(ut1), 1, n, length(ut1));
    B = zeros(length(ut1));
    for j = 1:length(t1)
        indk = find(t1(j) == ut1);
        indh = find(t2(j) == ut1);
        B(indk, indh) = -1; B(indh, indk) = -1;
    end
end


