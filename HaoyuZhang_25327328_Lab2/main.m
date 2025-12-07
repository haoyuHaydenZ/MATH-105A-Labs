% Q3
f  = @(x) x.^3/3 - x - 1/3;
p  = -0.347296;
I  = linspace(-0.5, 0.5, 1000);

figure; plot(I, f(I), 'LineWidth', 1.5); 
grid on; hold on;
yline(0,'k--');                          
xline(p,'r:');                           
plot(p, f(p), 'ro', 'MarkerSize', 6, 'MarkerFaceColor','r');
xlabel('x');
ylabel('f(x)'); 
title('f(x) = x^3/3 - x - 1/3 on [-0.5,0.5]');
legend('f(x)','y=0','x=p','location','best');

val_at_p = f(p);                %Verify the root numerically
fprintf('f(p) = %.10g\n', val_at_p);

% Q4 Bisection with error plot
f = @(x) x.^3/3 - x - 1/3;
a0 = -0.5; b0 = 0.5;         
tolF = 1e-4;                  
maxit = 100;

a = a0; b = b0; fa = f(a); fb = f(b);
if fa*fb > 0, error('f(a0) and f(b0) must have opposite signs.'); 
end

p_hist = [];
for k = 1:maxit
    p = (a + b)/2;
    p_hist(end+1,1) = p;    
    fp = f(p);

    if abs(fp) <= tolF, break; 
    end

    if fa*fp < 0
        b = p; fb = fp;
    else
        a = p; fa = fp;
    end
end

p_true = -0.347296;                 

n = (1:numel(p_hist))';
abs_err = abs(p_hist - p_true);
bound   = (b0 - a0) ./ (2.^n);     

figure;
semilogy(n, abs_err, 'o-', 'LineWidth', 1.2); 
hold on;
semilogy(n, bound,   '--', 'LineWidth', 1.2);
grid on;
xlabel('Iteration n'); 
ylabel('Error (log scale)');
title('|p_n - p| and bound (|b-a|/2^n) for bisection');
legend('|p_n - p|','(|b-a|)/2^n','Location','northeast');

% Q5 Fixed-point iteration with improved bound
f  = @(x) x.^3/3 - x - 1/3;
g  = @(x) (x.^3 - 1)/3;           
a = -0.5; b = 0.5;
K = 0.25;                         
p0 = 0;                          
tolF = 1e-4;
maxit = 200;

p = p0;
p_hist = []; f_hist = [];
for k = 1:maxit
    p = g(p);
    p_hist(end+1,1) = p;         
    f_hist(end+1,1) = f(p);
    if abs(f_hist(end)) <= tolF
        break;
    end
end

iters = numel(p_hist);
p_star = p_hist(end);
disp(['Fixed-point converged in ', num2str(iters), ...
      ' iterations. p ≈ ', num2str(p_star), ...
      ', |f(p)| = ', num2str(abs(f_hist(end)))]);

p_true = -0.347296;
n = (1:iters)';

abs_err = abs(p_hist - p_true);         
bound   = ((K.^n) / (1 - K)) * abs(g(p0) - p0); 

figure;
semilogy(n, abs_err, 'o-', 'LineWidth', 1.2); 
hold on;
semilogy(n, bound, '--', 'LineWidth', 1.2);
grid on;
xlabel('Iteration n');
ylabel('Error (log scale)');
title('|p_n - p| and theoretical bound (K^n/(1-K))|p_1 - p_0|');
legend('|p_n - p|','Bound','Location','northeast');

% Q6 Fixed-point with different initializations
f = @(x) x.^3/3 - x - 1/3;
g = @(x) (x.^3 - 1)/3;          
tolF = 1e-4;  maxit = 200;

p_true = -0.347296;              
p0s = linspace(-0.5, 0.5, 21);   

iters = zeros(size(p0s));
rate_est = NaN(size(p0s));

for i = 1:numel(p0s)
    p = p0s(i);
    p_hist = p; fh = f(p);
    for k = 1:maxit
        p = g(p);
        p_hist(end+1,1) = p;   
        fh = f(p);
        if abs(fh) <= tolF, break; 
        end
    end
    iters(i) = numel(p_hist)-1;

    if numel(p_hist) >= 6
        e = abs(p_hist - p_true);
        r = e(2:end)./e(1:end-1);           % e_{n+1}/e_n
        rate_est(i) = median(r(end-5:end-1));
    end
end

figure;
plot(p0s, iters, 'o-', 'LineWidth', 1.2); 
grid on;
xlabel('initialization p_0'); 
ylabel('iterations to |f(p_n)| \le 1e-4');
title('Fixed-point iterations vs initialization');

figure; 
plot(p0s, rate_est, 's-', 'LineWidth', 1.2); 
grid on;
xlabel('initialization p_0'); 
ylabel('estimated linear rate');
title('Estimated convergence rate vs initialization');
yline(0.25,'k--','K = sup|g''(x)| on [-0.5,0.5]');
yline((0.347296)^2,'r:','|g''(p)| \approx 0.1206');
