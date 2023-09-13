function [V, H, params, lucky] = lanczos_krylov(varargin)
%LANCZOS_KRYLOV  Krylov projection of a Hermitian matrix A.
%
% [V, H] = LANCZOS_KRYLOV(A, B) construct the Krylov
%     subspace spanned by [B, A*B]. The matrix V is an orthogonal
%     basis for this space, H is block tridiagonal and satisfies
%
%        A * V = V * H                                          (1)
%
% [V, H, params] = LANCZOS_KRYLOV(V, H, PARAMS) enlarges the
%     Krylov subspace generated with a previous call to LANCZOS_KRYLOV. 
%     
% Note: lucky is a flag parameter that highlight the lucky breakdown of lanczos


if nargin ~= 2 && nargin ~= 3
    error('Called with the wrong number of arguments');
end

if nargin == 2
    % Start to construct the Krylov space
    [V, H, params, lucky] = lanczos_krylov_start(varargin{:});
else
    % Enlarge the space that was previously built
    [V, H, params, lucky] = lanczos_krylov_extend(varargin{:});	
end

end

function [V, H, params, lucky] = lanczos_krylov_start(A, b)

if ~isstruct(A)
    m = size(A, 1);
    n = size(A, 2);
    
    if m ~= n
        error('The matrix A should be square');
    end
    
    if n ~= size(b, 1)
        error('The block vector b has wrong number of rows');
    end
end

bs = size(b, 2);

% Construct a basis for the column span of b
[V, ~] = qr(b, 0);

H = zeros(bs, 0);

[V, H, w, lucky] = add_inf_pole (V, H, A, V);

% Save parameters for the next call
params = struct();
params.last = w;
params.A = A;
end

function [V, H, params, lucky] = lanczos_krylov_extend(V, H, params)
w  = params.last;
A  = params.A;

[V, H, w, lucky] = add_inf_pole (V, H, A, w);

params.last = w;
end

%
% Utility routine that adds an infinity pole to the space. The vector w is
% the continuation vector.
%
function [V, H, w, lucky] = add_inf_pole(V, H, A, w)
lucky_tol = 1e-8; % tolerance for detecting lucky breakdowns TODO handle more general breakdowns
lucky = false;
bs = size(w, 2);

if isstruct(A)
    w = A.multiply(1.0, 0.0, w);
else
    w = A * w;
end

% Enlarge H and K
H(size(H, 1) + bs, size(H, 2) + bs) = 0;

% Perform orthogonalization with modified Gram-Schimidt
[w, H(max(1,end - 3*bs +1):end-bs, end-bs+1:end)] = mgs_orthogonalize(V, w);

[w, H(end-bs+1:end, end-bs+1:end)] = qr(w, 0);
if norm(H(end-bs+1:end, end-bs+1:end), 'fro') < lucky_tol
	lucky = true;
end
if size(V, 2) == bs
	V = [V, w];
else
	V(:, 1:bs) = w;
	V = V(:, [bs+1:2 * bs, 1:bs]);
end

end

%
% Modified Gram-Schmidt orthogonalization procedure.
%
% Suggested improvements: work with block-size matrix vector products to
% get BLAS3 speeds.
%
function [w, h] = mgs_orthogonalize(V, w)
    h = V' * w;
    w = w - V * h;
    h1 = V' * w;
    h = h + h1;
    w = w - V * h1;
end

