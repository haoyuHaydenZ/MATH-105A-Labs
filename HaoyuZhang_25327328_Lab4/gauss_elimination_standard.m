function x = gauss_elimination_standard(A, b)
% Inputs:  A (n-by-n), b (n-by-1)
% Output:  x (n-by-1)

    [n,m] = size(A);
    if n ~= m, error('A must be square'); 
    end
    if ~isvector(b) || length(b) ~= n, error('b must be n-by-1'); 
    end
    b = b(:);

    U = double(A);
    y = double(b);

    tol = 1e-15;
    scale = norm(U, inf);

    for k = 1:n-1
        pivot = U(k,k);
        if abs(pivot) <= tol * max(1, scale)
            error('Zero/near-zero pivot at step %d (no pivoting).', k);
        end
        for i = k+1:n
            m_ik = U(i,k)/pivot;
            if m_ik ~= 0
                U(i, k:n) = U(i, k:n) - m_ik * U(k, k:n);
                y(i)      = y(i)      - m_ik * y(k);
            end
        end
    end

    if abs(U(n,n)) <= tol * max(1, scale)
        error('Zero/near-zero final pivot.');
    end

    x = zeros(n,1);
    for i = n:-1:1
        s = U(i, i+1:end) * x(i+1:end);
        x(i) = (y(i) - s) / U(i,i);
    end
end
