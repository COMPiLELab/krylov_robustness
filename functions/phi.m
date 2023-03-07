function nrm = phi(a, A, X, S, r)
    fAX = expm(full(A + X)); v = diag(fAX);
    nrm = norm(r(S) - a .* v(S), 2)^2;
end
