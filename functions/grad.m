function [g,G] = grad(a, A, X, S, EdgesPattern, r)
    fAX = expm(full(A + X)); v = diag(fAX);
    g = 2 * (a * v(S) - r(S))' * v(S);
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
            
            G(i,j) = 2 * a * deB(S)' * (a * v(S) - r(S));
        end
        
    end
G = G + G';
end
