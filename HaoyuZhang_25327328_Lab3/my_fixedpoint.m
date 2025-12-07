function [p, fp, iter] = my_fixedpoint(g, p0, tol, max_iter)

    if nargin < 3
        tol = 1e-6;
    end
    if nargin < 4
        max_iter = 100;
    end

    p = p0;

    for k = 1:max_iter
        p_next = g(p);
        
        if abs(p_next - p) < tol
            p = p_next;
            fp = g(p) - p;
            iter = k;
            fprintf('Converged in %d iterations.\n', k);
            return
        end
        
        p = p_next;
    end

    fp = g(p) - p;
    iter = max_iter;
    fprintf('Max iterations reached. Approx fixed point = %.8f\n', p);
end
% %->start of a format specifier; .8->show 8 digits after the decimal
% point; f->display as a floating-point number (decimal format); \n ->new
% line