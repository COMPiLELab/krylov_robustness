datasets = {'london', 'DC_street_graph'};

for i =  1:length(datasets)
	filename = sprintf('../Datasets/%s.mat', datasets{i});
	load(filename) 

    A = spones(Problem.A+Problem.A'); 
	A = A - spdiags(diag(A), 0, size(A, 1), size(A, 2));

    %%% get max connected component
    ind = max_connected_component(A);
    A = A(ind,ind);
    n = size(A,1);
    if n < 5000
        centrality = compute_centrality(A,'exp');
    else
        centrality = compute_centrality(A,'eig');
    end
	
	for j = 1:n
		E1 = find_top_missing_edges(A, centrality, j);
		oind = find(E1(:, 1) >= E1(:, 2));
		E1(oind, :) = [E1(oind, 2), E1(oind, 1)]; 
		E2 = brute_force(A, centrality, j);
		if norm(E1-E2) > 0
			disp('Discrepancy between the two set of edges\n')
			keyboard
		end
	end
	pause
end


function E = brute_force(A, centrality, num)
	centrality = centrality(:);
	S = centrality * centrality';
	S = triu(S);
 	[ ~, ind ] = sort(S(:), 1, 'descend');
     [ I, J ] = ind2sub(size(S), ind);
     E = [];
     found_edges = 0;
     i = 1;
     while found_edges < num && i <= length(I)
         if A(I(i),J(i)) == 0 && I(i)~=J(i) 
             E = [E; [I(i) J(i)]];
             found_edges = found_edges + 1;
         end
         i = i+1;
        % if mod(i,10) == 0
        %     fprintf('i = %d\n',i)
        % end
     end
end

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
