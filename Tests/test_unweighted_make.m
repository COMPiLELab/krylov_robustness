% Test the performances of the various greeedy schemes for the break problem on either benchmark or transportation graphs (determined by datasets{1} or dataset{2} in the for loop)
% with a fixed budget k=50 of edges to be added
clear all
addpath ../functions
addpath ../'MIOBI Codes'/

debug = false;

datasets{1} = dir('../datasets_paper/Misc/*');
datasets{2} = dir('../datasets_paper/Transport/*');
column_names = {'method', 'dataset', 'n', 'm', ...
                'searchspace_size', 'centrality_order', 'time', 'tr_variation','budget_size'};
Results_TAB = table([],[],[],[],[],[],[],[],[],'VariableNames',column_names);
 
tol = 1e-6; 
poles = inf; 
it = 100; 
debug = 0;
k = 50; % budget size      
Q = 250; % parameter that tunes the size of the search space in greedy_krylov_make	  
intersection_sizes = zeros(0, 4); % store the cardinality of the edges that are in common among the various methods (3 pairs and 1 triple)

for ii = 1:2 % Loop over benchmarks: miscellaneous graphs and transportation networks

	data = datasets{ii}; 
	l = 3;
	
	for i =  l:length(data) %start from 3 as first two are '.' and '..'
		if ii == 1
			filename = sprintf('../datasets_paper/Misc/%s', data(i).name);
		else	
			filename = sprintf('../datasets_paper/Transport/%s', data(i).name);			
		end	
		name = split(data(i).name,'.');
		name = string(name{1});
		load(filename); 
		A = Problem.A;

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
		tic;
		centrality = compute_centrality(A,'eig');
		centrality_order = {'mult','min'}; % first and second orders based on the eigenvector centrality (see the manuscript)
		time_centrality = toc;
		
		fprintf('Dataset: %s\t n: %d\t budget: %d\t ||exp(A)|| = %.1e\n', name, n, k, nrm);

        %% GREEDY_KRYLOV_MAKE (centrality order employed = 2)            
        method  = "GREEDY_KRYLOV_MAKE";
        Q = min(nnz(A)/2 - k, Q);
        tic;
        size_search_space = Q + k;
        [edges, delta_trace, A_new] = greedy_krylov(A, k, Q, centrality, centrality_order{2}, tol * nrm, it, poles, debug,'make');
        time_greedy_krylov_make = toc;
        time_greedy_krylov_make = time_greedy_krylov_make + time_centrality;
        Results_TAB = [Results_TAB; ...
        {method, name,n,nnz(A)/2,size_search_space,...
        string(centrality_order{2}),time_greedy_krylov_make,delta_trace/trexp,k}];
    
		if debug 
			edadd = nnz(spones(spones(A) - spones(A_new)));
			fprintf('Number of edges added by GREEDY_KRYLOV_MAKE = %d\n', edadd/2);
		end

	    
	    %% MIOBI 25 eigs
	    method = "MIOBI";
	    num_eig_miobi = 25; % num of eigenvectors to use in MIOBI
	    tic;
		[~, ~, Aorg] = MIOBIMakeEdge(A,k,num_eig_miobi);
		[t1, t2, ~] = find(spones(spones(A) - spones(Aorg)));
		[U, B] = edge2low_rank([t1, t2], n);
	    delta_trace = trace_fun_update(A, full(U), B, tol * nrm);
		time_miobi = toc;
		time_miobi = time_miobi + time_centrality;
	    
	    [~,ind] = sort(centrality,'descend');
	    search_space_miobi = (dmax^2 - nnz(A(ind(1:dmax),ind(1:dmax))) - dmax)/2;
	    
	    Results_TAB = [Results_TAB; ...
	        {method, name,n,nnz(A)/2,search_space_miobi,...
	        "--",time_miobi,delta_trace/trexp,k}];

		if debug 
			fprintf('MIOBI done\n')
		end

	    %% EIGENV (Arrigo&Benzi, SISC 2016)
	    method = "EIGENV";
	    tic; 
	    [~,ind] = sort(centrality, 'descend');
	    ind_top_nodes = ind(1:ceil(n/5));
	    Asmall = A(ind_top_nodes,ind_top_nodes); 
	    if size(Asmall,1)^2-nnz(Asmall)-size(Asmall,1)<2*k
	        Asmall = A;
	        ind_top_nodes = [1:n]';
	    end
	    
	    EE = find_top_missing_edges(Asmall,centrality(ind_top_nodes),k,'mult');
	    HURISTICEdges = [ind_top_nodes(EE(:,1)) ind_top_nodes(EE(:,2))];
	    
	    [U, B] = edge2low_rank(HURISTICEdges, n);
		delta_trace = trace_fun_update(A, full(U), B, tol * nrm);
		time_eigenv = toc;
		time_eigenv = time_eigenv + time_centrality;
		Results_TAB = [Results_TAB; ...
	        {method, name,n,nnz(A)/2,k,"mult",time_eigenv,delta_trace/trexp,k}];
		if debug 
			fprintf('EIGENV done\n')
		end
	   
		% Compute numbers of common edges
		edges = sort(edges, 2);
		edges_miobi = sort([t1, t2], 2);
		HURISTICEdges = sort(HURISTICEdges, 2);
		inter = zeros(1, 4);
		tmp = intersect(edges, edges_miobi, 'rows'); % intersection between GREEDY_KRYLOV_BREAK and MIOBI
		inter(1) = size(tmp, 1);
		tmp = sort(tmp, 2);
		tmp = intersect(tmp, HURISTICEdges, 'rows'); % further intersection with EIGENV
		inter(4) = size(tmp, 1);
		tmp = intersect(edges, HURISTICEdges, 'rows'); % intersection between GREEDY_KRYLOV_BREAK and EIGENV	
		inter(2) = size(tmp, 1);	
		tmp = intersect(edges_miobi, HURISTICEdges, 'rows'); % intersection between MIOBI and EIGENV
		inter(3) = size(tmp, 1);				
		intersection_sizes = [intersection_sizes; inter];
	   
	    filepath = sprintf('../Results/results_unweighted_make_%s.csv', string(date));
	    writetable(Results_TAB,filepath)
		disp(Results_TAB(Results_TAB.dataset == name,:))
	   	fprintf('Number common edges: GKB & MIOBI = %d \t GKB & EIGENV = %d, \t MIOBI & EIGENV = %d \t ALL = %d\n', intersection_sizes(end, :));
		    
	end
end

dlmwrite('../Results/results_intersection_make.dat', intersection_sizes, '\t');


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
%---------------------------------------------------------
function [U, B] = edge2low_rank(E, n)
% Build the low-rank factorization that adds the edges in E
 	t1 = E(:,1);
    t2 = E(:,2);
    ut1  = unique([t1; t2]);
    U = sparse(ut1, 1:length(ut1), 1, n, length(ut1));
    B = zeros(length(ut1));
    for j = 1:length(t1)
        indk = find(t1(j) == ut1);
        indh = find(t2(j) == ut1);
        B(indk, indh) = 1; B(indh, indk) = 1;
    end
end




