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
