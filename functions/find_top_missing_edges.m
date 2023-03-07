function E = find_top_missing_edges(A, centrality, num, order) 
% Returns the indices of the top 'num' missing edges
%---------------------------INPUT-------------------------------------------------------------------------------------
% 
% A				adjacency matrix of the graph
% centrality	centrality measure associated with the nodes
% num			number of edges with top centrality measure to return
% order			(optional) ordering strategy on the set of edges: 
%				'mult' --> multiplication of nodes centrality
%				'min'  --> first compare the minimum, if equal compare the maximum of the node centralities (default)
%
%---------------------------------------------------------------------------------------------------------------------
%
	if ~exist('order', 'var')
		order = 'min';
	end
	centrality = centrality(:);	
	[sc, indC] = sort(centrality, 'descend');
	switch order
		case 'mult' % Ordering based on multiplication of node centralities
			if (size(A, 1)^2 - nnz(A) - size(A, 1))/2 <= num
				S = centrality(indC) * centrality(indC)';
				S = triu(S);
				[ ~, ind ] = sort(S(:), 1, 'descend');
				[ I, J ] = ind2sub(size(S), ind); 
				I = indC(I); J = indC(J);
				return
			end
			% Compute a minimum value for N (there has to be enough missing edges)
			min_N = 2;
			len = 0;
			while len < num
				len = len + sum( A( indC(1 : min_N - 1), indC(min_N) ) == 0 );
			 	min_N = min_N + 1;
			end
			min_N = min_N - 1;

			 N = sum(sc(1) * sc > sc(min_N)^2);
			 S = centrality(indC(1:N)) * centrality(indC(1:N))';
			 S = triu(S);
			 [ ~, ind ] = sort(S(:), 1, 'descend');
			 [ I, J ] = ind2sub(size(S), ind); 
			 I = indC(I); J = indC(J);
			
			 E = [];
			 found_edges = 0;
			 i = 1;
			 while found_edges < num && i <= length(I)
				 if A(I(i),J(i)) == 0 && I(i)~=J(i) 
				     E = [E; [I(i) J(i)]];
				     found_edges = found_edges + 1;
				 end
				 i = i + 1;
			 end
		case 'min' % Ordering based on minum and maximum of node centralities
			E = [];
			j = 2;
			while size(E, 1) < num
				ind = find(A( indC(1 : j - 1), indC(j) ) == 0);
				E = [E; [indC(ind)], indC(j) * ones(length(ind), 1)];
				j = j + 1;
			end
			if size(E, 1) > num
				E = E(1:num, :);
			end
	end        
end
