%% Problem 2 Setup:
f  = @(x) x - 0.8 - 0.2*sin(x);
df = @(x) 1 - 0.2*cos(x);
g  = @(x) 0.8 + 0.2*sin(x);   % |g'(x)| = 0.2*cos(x) <= 0.2 < 1 => contraction

tol = 1e-4;

% (a)
[a, b] = deal(0, pi/2);
[root_bis, froot_bis, it_bis] = my_bisection(f, a, b, tol);
fprintf('(a) Bisection: root ≈ %.7f,  f(root)=%.2e,  iters=%d\n', root_bis, froot_bis, it_bis);

% (b)
p0 = pi/4;
[p_fp, fp_resid, it_fp] = my_fixedpoint(g, p0, tol);
fprintf('(b) Fixed point: root ≈ %.7f,  residual g(p)-p=%.2e,  iters=%d\n', p_fp, fp_resid, it_fp);

% (c)
root_newton = my_newtons(f, df, p0, tol, 50);
fprintf('(c) Newton: root ≈ %.7f,  f(root)=%.2e\n', root_newton, f(root_newton));

%% Problem 3 Setup:
f  = @(x) x - 0.8 - 0.2*sin(x);
df = @(x) 1 - 0.2*cos(x);
g  = @(x) 0.8 + 0.2*sin(x);
tol = 1e-4;

% Bisection:
a = 0; b = pi/2;
fa = f(a); fb = f(b);
p_bis = []; f_bis = [];
for k = 1:100
    c = (a + b)/2;
    fc = f(c);
    p_bis(end+1) = c;
    f_bis(end+1) = fc;
    if abs(fc) < tol || (b - a)/2 < tol, break; 
    end
    if fa*fc < 0, b=c; fb=fc; 
    else, a=c; fa=fc; 
    end
end

% Fixed Point:
p = pi/4;
p_fp = p; f_fp = f(p);
for k = 1:100
    p_next = g(p);
    p_fp(end+1) = p_next;
    f_fp(end+1) = f(p_next);
    if abs(p_next - p) < tol, break;
    end
    p = p_next;
end

% Newton's
p = pi/4;
p_newt = p; f_newt = f(p);
for k = 1:50
    p_next = p - f(p)/df(p);
    p_newt(end+1) = p_next;
    f_newt(end+1) = f(p_next);
    if abs(p_next - p) < tol, break; 
    end
    p = p_next;
end

% Plotting:
figure; 
hold on;
semilogy(0:length(f_bis)-1, abs(f_bis), 'o-','DisplayName','Bisection');
semilogy(0:length(f_fp)-1, abs(f_fp), 's-','DisplayName','Fixed Point');
semilogy(0:length(f_newt)-1, abs(f_newt), '^-','DisplayName','Newtons');
xlabel('Iteration n'); 
ylabel('|f(p_n)|');
legend('Location','southwest'); 
grid on;
title('|f(p_n)| per iteration with three methods');

%% Problem 4 Setup: 
f  = @(x) x.^2 - 10*cos(x);
df = @(x) 2*x + 10*sin(x);
tol = 1e-5;  max_iter = 100;

p0s = [-26 -25 -24 24 25 26];
for j = 1:numel(p0s)
    p = p0s(j);
    for k = 1:max_iter
        p_new = p - f(p)/df(p);
        if abs(p_new - p) < tol, break; 
        end
        p = p_new;
    end
    fprintf('p0=%3d  ->  root ≈ %.9f   iters=%2d   |f(root)|=%.2e\n', ...
            p0s(j), p_new, k, abs(f(p_new)));
end
