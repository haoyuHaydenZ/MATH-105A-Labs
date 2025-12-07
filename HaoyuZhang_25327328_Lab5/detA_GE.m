function d = detA_GE(A)
%   d = detA_GE(A) computes det(A) for a real square matrix using
%   row-reduction to upper-triangular form. Only row swaps are used;
%   elimination adds multiples of rows (which does not change det).
%   Row swap flips the sign of the determinant.
%   If a pivoted column is entirely ~0, the matrix is singular -> det = 0.

%   Example:
%       A = [2 -1 0; -1 2 -1; 0 -1 2];
%       detA_GE(A)

    if ndims(A) ~= 2 || size(A,1) ~= size(A,2)
        error('detA_GE:InputNotSquare', 'Input must be a square matrix.');
    end
    n = size(A,1);

    if n == 0
        d = 1;
        return
    elseif n == 1
        d = A(1,1);
        return
    end

    U = A;          % upper triangular
    sgn = 1;        % sign from row swaps
    tol = eps * max(1, norm(A,'fro'));

    for k = 1:n-1
        % pivot: index of largest |entry| in column k, rows k:n
        [~, pRel] = max(abs(U(k:n,k)));
        p = pRel + k - 1;

        if abs(U(p,k)) <= tol
            d = 0;
            return
        end

        % row swap if needed (flip sign)
        if p ~= k
            tmp     = U(k,:); U(k,:) = U(p,:); U(p,:) = tmp;
            sgn     = -sgn;
        end

        % eliminate below pivot
        for i = k+1:n
            m = U(i,k) / U(k,k);
            if m ~= 0
                % row_i := row_i - m * row_k; That does not change det.
                U(i, k:n) = U(i, k:n) - m * U(k, k:n);
            end
        end
    end

    if abs(U(n,n)) <= tol
        d = 0;
        return
    end

    % determinant = sign * product of diagonal of U
    d = sgn * prod(diag(U));
end
