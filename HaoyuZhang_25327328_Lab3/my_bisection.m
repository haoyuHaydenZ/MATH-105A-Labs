function [root, froot, iter] = my_bisection(f, a, b, tol, max_iter)

    if nargin < 4
        tol = 1e-6;
    end
    if nargin < 5
        max_iter = 100;
    end

    fa = f(a);
    fb = f(b);

    if fa * fb > 0
        error('f(a) and f(b) must have opposite signs.');
    end

    for k = 1:max_iter
        c = (a + b)/2;
        fc = f(c);

        if abs(fc) < tol || (b - a)/2 < tol
            root = c;
            froot = fc;
            iter = k;
            fprintf('Converged in %d iterations.\n', k);
            return
        end

        if fa * fc < 0
            b = c;
            fb = fc;
        else
            a = c;
            fa = fc;
        end
    end

    root = (a + b)/2;
    froot = f(root);
    iter = max_iter;
    fprintf('Max iterations reached and the approx root = %.8f\n', root);
end
