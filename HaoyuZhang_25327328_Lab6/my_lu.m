function [L, U] = my_lu(A)
%   performs an LU factorization of a square matrix A such that A = L*U,
%   L: lower-triangular; U: upper-triangular
%       A : n-by-n matrix
%   OUTPUT - L : n-by-n unit lower-triangular matrix; U : n-by-n upper-triangular matrix
%   If a zero pivot occurs, the routine stops and displays an error message.
%   Reference: (lu pseudo.pdf from the Demos and Notes folder) Algorithm 6.4, Numerical Analysis (9th Edition), Burden & Faires
[n, m] = size(A);
if n ~= m
    error('Matrix A must be square.');
end

L = eye(n);
U = zeros(n);

% Step 1: first pivot
if A(1,1) == 0
    error('Factorization impossible: zero pivot at (1,1).');
end
L(1,1) = 1;
U(1,1) = A(1,1);

% Step 2: first row of U and first column of L
for j = 2:n
    U(1,j) = A(1,j) / L(1,1);   
    L(j,1) = A(j,1) / U(1,1);
end

% Step 3–5: inner loops
for i = 2:n-1
    % Step 4: compute U(i,i)
    sumLU = 0;
    for k = 1:i-1
        sumLU = sumLU + L(i,k)*U(k,i);
    end
    U(i,i) = A(i,i) - sumLU;
    L(i,i) = 1;

    if U(i,i) == 0
        error('Factorization impossible: zero pivot at (%d,%d).', i, i);
    end

    % Step 5: compute remaining entries in row i of U and column i of L
    for j = i+1:n
        sum1 = 0;
        for k = 1:i-1
            sum1 = sum1 + L(i,k)*U(k,j);
        end
        U(i,j) = (A(i,j) - sum1) / L(i,i);

        sum2 = 0;
        for k = 1:i-1
            sum2 = sum2 + L(j,k)*U(k,i);
        end
        L(j,i) = (A(j,i) - sum2) / U(i,i);
    end
end

% Step 6: final pivot (n,n)
sumEnd = 0;
for k = 1:n-1
    sumEnd = sumEnd + L(n,k)*U(k,n);
end
U(n,n) = A(n,n) - sumEnd;
L(n,n) = 1;

if U(n,n) == 0
    warning('A = LU but A is singular (U(n,n)=0).');
end

end
