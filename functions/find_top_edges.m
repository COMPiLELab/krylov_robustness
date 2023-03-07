function E = find_top_edges(A, centrality, num, order) 
% Returns the indices of the top 'num' existing edges
%---------------------------INPUT---------------------------------------------------------------------------
% 
% A				adjacency matrix of the graph
% centrality	centrality measure associated with the nodes
% num			number of edges with top centrality measure to return
% order			(optional) ordering strategy on the set of edges: 
%				'mult' --> multiplication of nodes centrality (default)
%				'min'  --> first compare the minimum, if equal compare the maximum of the node centralities
%
%-----------------------------------------------------------------------------------------------------------
%
	if ~exist('order', 'var')
		order = 'mult';
	end
	[I, J] = find(tril(A, -1));
	E = [I, J];
	if length(I) < num
		warning('FIND_TOP_EDGES:: there are not enough edges in the graph');
	end
	switch order
		case 'mult'
			c = centrality(I) .* centrality(J);
			[~, ind] = sort(c, 'descend');
			E = E(ind(1:num), :);       
		case 'min'
			sc = sort(centrality, 'descend');
			scores = zeros(size(E, 1),1);
			for h = 1:size(scores, 1)
				c1 = find(sc == centrality(I(h)), 1);
				c2 = find(sc == centrality(J(h)), 1);
				mn = min(c1, c2);
				mx = max(c1, c2);
				scores(h) = mx * (mx - 1) / 2 + mn;  
			end
			[~, ind] = sort(scores, 'ascend');
			E = E(ind(1:num), :);       
	end 
end
