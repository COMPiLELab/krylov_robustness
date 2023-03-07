function [g,G] = grad_log(a, A, X, difAX, S, EdgesPattern, r, mu)
%---------------------------------------------------------------------------
% Evaluate the gradient of the objective function with log-barrier
%
%   a, X		parameters of the modification
%   A 			adjacency matrix of the original graph
%   difAX               starting (centrality) score
%   S			set of nodes where to evaluate the fitting of scores
%   EdgesPattern 	sparsity pattern of the perturbation
%   r			target (centrality) score
%   mu		 	multiplicative parameter of the log-barrier
%---------------------------------------------------------------------------
edg = sub2ind(size(A),EdgesPattern(:,1),EdgesPattern(:,2));

%fAX = expm(full(A + X)); v = diag(fAX);
g = 2 * (a * difAX(S) - r(S))' * difAX(S);
G = zeros(size(A));
n = length(r);



for k = 1: size(EdgesPattern, 1)
    i = EdgesPattern(k,1);
    j = EdgesPattern(k,2);
    
    if i < j % works only for symmetric adjacency matrices
        Eij = zeros(size(A)); Eij(i,j) = 1;
        B = [A+X Eij; zeros(size(A)) A'+X'];
        eB = expm(full(B));
        deB = diag(eB(1:n,n+1:end));
        
        G(i,j) = 2 * a * deB(S)' * (a * difAX(S) - r(S));
    end
end
    
G = G + G';

g = g - mu./a;

temp = zeros(size(A)); 
temp(edg) = 1./(A(edg)+X(edg));

G = G - mu .* temp;

end
