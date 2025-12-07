clear; clc;

A = [ 2  1  1;
      4 -6  0;
     -2  7  2 ];

B = [1e-4 1;
     1     1];

n = 30;
C = tril(-ones(n), -1) + eye(n);   
C(:, end) = 1;                   

Ms    = {A, B, C};
names = {'A','B','C'};

gf_std = zeros(1,numel(Ms));
gf_pp  = zeros(1,numel(Ms));
for t = 1:numel(Ms)
    gf_std(t) = growth_factor(Ms{t}, 'standard');
    gf_pp(t)  = growth_factor(Ms{t}, 'partial');
end

fprintf('\nTable 1: Growth Factors for solving Mx = 0\n');
fprintf('--------------------------------------------------------\n');
fprintf('%-16s %-14s %-14s %-14s\n', '', 'A', 'B', 'C');
fprintf('--------------------------------------------------------\n');
fprintf('%-16s %-14.6e %-14.6e %-14.6e\n', 'Standard Pivoting', gf_std(1), gf_std(2), gf_std(3));
fprintf('%-16s %-14.6e %-14.6e %-14.6e\n', 'Partial Pivoting',  gf_pp(1),  gf_pp(2),  gf_pp(3));
fprintf('--------------------------------------------------------\n\n');

% Is the output of GE the solution to Mx = 0?
tol = 1e-12;
for t = 1:numel(Ms)
    M = Ms{t};
    b = zeros(size(M,1),1);
    try
        x_std = gauss_elimination_standard(M, b);
        x_pp  = gauss_elimination_partial(M, b);
    catch ME
        warning('Solver error for %s: %s', names{t}, ME.message);
        continue
    end
    fprintf('Matrix %s: ||x_std||=%.2e, ||x_pp||=%.2e, det≈%.2e\n', ...
        names{t}, norm(x_std), norm(x_pp), det(M));
    if norm(x_std) < tol && norm(x_pp) < tol
        fprintf('  -> Output is (numerically) the zero vector: YES\n');
    else
        fprintf('  -> Output deviates from zero due to rounding.\n');
    end
end

function rho = growth_factor(M, method)
    U = double(M);
    [n,~] = size(U);
    tol = 1e-15;
    scale = norm(U,inf);

    max_init   = max(abs(U(:)));
    max_during = max_init;

    for k = 1:n-1
        if strcmpi(method,'partial')
            [~, idx] = max(abs(U(k:n, k)));
            p = k + idx - 1;
            if p ~= k
                U([k p], :) = U([p k], :);
            end
        elseif ~strcmpi(method,'standard')
            error('Unknown method: %s', method);
        end

        pivot = U(k,k);
        if abs(pivot) <= tol * max(1, scale)
            rho = max_during / max(max_init, eps);
            return;
        end
        for i = k+1:n
            m = U(i,k)/pivot;
            if m~=0
                U(i,k:n) = U(i,k:n) - m*U(k,k:n);
            end
        end
        max_during = max(max_during, max(abs(U(:))));
    end
    rho = max_during / max(max_init, eps);
end
