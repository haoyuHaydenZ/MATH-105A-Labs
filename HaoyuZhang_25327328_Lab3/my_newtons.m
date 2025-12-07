function root = my_newtons(f, df, p0, tol, max_iter)

for k = 1:max_iter
    p1 = p0 - f(p0) / df(p0);    
    if abs(p1 - p0) < tol      
        root = p1;
        fprintf('Converged after %d iterations.\n', k);
        return;
    end
    p0 = p1;
end

warning('Newton''s Method did not converge after %d iterations.', max_iter);
root = p1;
end
