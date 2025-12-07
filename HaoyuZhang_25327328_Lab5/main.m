function main_det_experiments()
% MAIN_DET_EXPERIMENTS
%  (1) Timing: n = 1..10, A = randn(n,n); compare detA_recursive vs detA_GE
%  (2) Round-off sensitivity: n = 1000, D = (1/100)*I; compute det and discuss invertibility

    close all; clc;

    %% (1) TIMING
    rng(42);                    
    ns   = 1:10;
    tRec = nan(size(ns));
    tGE  = nan(size(ns));
    dRec = nan(size(ns));
    dGE  = nan(size(ns));

    fprintf('--- Timing experiment (randn matrices) ---\n');
    for idx = 1:numel(ns)
        n = ns(idx);
        A = randn(n,n);

        % recursive timing
        tic;
        dRec(idx) = detA_recursive(A);
        tRec(idx) = toc;

        % GE timing
        tic;
        dGE(idx)  = detA_GE(A);
        tGE(idx)  = toc;

        fprintf('n=%2d | time(rec)=%.4fs  time(GE)=%.4fs  | det diff=%.3e\n', ...
            n, tRec(idx), tGE(idx), abs(dRec(idx)-dGE(idx)));
    end

    % plot: time vs n
    figure('Color','w');
    plot(ns, tRec, '-o', 'LineWidth',1.25, 'MarkerSize',6); hold on;
    plot(ns, tGE , '-s', 'LineWidth',1.25, 'MarkerSize',6);
    ylabel('Time (seconds)'); xlabel('n');
    title('Determinant timing: recursive vs Gaussian elimination');
    legend({'detA\_recursive','detA\_GE'}, 'Location','northwest');
    grid on;
    set(gca,'YScale','log'); %(makes growth difference obvious)
    saveas(gcf, 'timing_recursive_vs_GE.png');

    %% (2) ROUND-OFF SENSITIVITY
    fprintf('\n--- Round-off sensitivity (diagonal D) ---\n');
    n = 1000;
    lambda = 1/100;                 
    D = diag(lambda*ones(n,1));

    % GE-based determinant (will underflow in double precision)
    dD_GE = detA_GE(D);

    % True value in exact arithmetic: lambda^n
    true_det = lambda^n;              

    % det(D) = prod(diag) = lambda^n, so log|det(D)| = n*log(lambda), sign = +1
    logAbsDet = n*log(lambda);
    signDet   = 1;

    fprintf('n = %d, lambda = %.2g\n', n, lambda);
    fprintf('detA_GE(D)        = %.5g   (double precision; underflow expected)\n', dD_GE);
    fprintf('sign(det(D))      = %d\n', signDet);
    fprintf('log(|det(D)|)     = %.5f  ->  |det(D)| = exp(logAbsDet) ~ 10^{%.1f}\n', ...
            logAbsDet, logAbsDet/log(10));

    % D is diagonal with nonzero diagonal entries => mathematically invertible.
    if lambda ~= 0
        hasInverse = true;
    else
        hasInverse = false;
    end
    fprintf('Is D invertible?  %s (yes, because all diagonal entries are nonzero)\n', ...
            ternary(hasInverse,'YES','NO'));
    fprintf('Note: The exact determinant is lambda^n = (1/100)^{1000} = 10^{-2000},\n');
    fprintf('      which is far below double-precision min (~5e-324), so it underflows to 0.\n');

end
