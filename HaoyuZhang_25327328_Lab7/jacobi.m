function [x, k] = jacobi(A, b, x0, tol, maxit)
%Solve A x = b using the Jacobi iterative method.
%[x, k] = jacobi(A, b, x0, tol, maxit)
%   A      : n x n coefficient matrix
%   b      : right-hand side vector
%   x0     : initial guess
%   tol    : relative tolerance for stopping (e.g. 1e-8)
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
        error('Jacobi:zeroDiagonal', 'Matrix has zero diagonal entries.');
    end

    x  = x0(:);                 %make sure the column vector
    x_new = x;

    D = diag(diag(A));
    R = A - D;                  %R = L + U

    for k = 1:maxit
        x_new = D \ (b - R * x);

        %Stopping criterion
        if norm(x_new - x, inf) <= tol * max(1, norm(x_new, inf))
            x = x_new;
            return
        end

        x = x_new;
    end

    %max iterations reached
    warning('Jacobi:maxit', 'Maximum number of iterations reached.');
end
