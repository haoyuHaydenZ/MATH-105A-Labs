function [x, k] = gauss_seidel(A, b, x0, tol, maxit)
%GAUSS_SEIDEL   Solve A x = b using the Gauss-Seidel method
%[x, k] = gauss_seidel(A, b, x0, tol, maxit)
%   A      : n x n coefficient matrix
%   b      : right-hand side vector
%   x0     : initial guess
%   tol    : relative tolerance (e.g. 1e-8)
%   maxit  : maximum number of iterations
%   x      : approximate solution
%   k      : number of iterations performed

    if nargin < 3 || isempty(x0)
        x0 = zeros(size(b));
    end
    if nargin < 4 || isempty(tol)
        tol = 1e-8;
    end
    if nargin < 5 || isempty(maxit)
        maxit = 1000;
    end

    if any(abs(diag(A)) < eps)
        error('GaussSeidel:zeroDiagonal', 'Matrix has zero diagonal entries.');
    end

    x = x0(:);                 %column vector
    n = length(b);

    for k = 1:maxit
        x_old = x;

        for i = 1:n
            %Use latest values x(1:i-1) and old values x_old(i+1:n)
            s1 = A(i,1:i-1) * x(1:i-1);
            s2 = A(i,i+1:n) * x_old(i+1:n);
            x(i) = (b(i) - s1 - s2) / A(i,i);
        end

        %the Stopping criterion
        if norm(x - x_old, inf) <= tol * max(1, norm(x, inf))
            return
        end
    end

    warning('GaussSeidel:maxit', 'Maximum number of iterations reached.');
end
