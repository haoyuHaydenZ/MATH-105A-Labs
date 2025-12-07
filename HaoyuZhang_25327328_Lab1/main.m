%% Main File
xx = 5/7;
yy = 1/3;
uu = 0.714215;
ww = 0.111111 * 10^-4;
kk = 5;

%% Problem 3: Compute values using 5 digit chopping
% (a) x ⊕ y
aa = xx + yy;
[ahat, ~] = fl_chop(aa, kk);

% (b) x ⊖ u
bb = xx - uu;
[bhat, ~] = fl_chop(bb, kk);

% (c) (x ⊖ u) ⊘ w  (DO chop after subtraction, then after division)
[tempDiff_chop, ~] = fl_chop(xx - uu, kk);      % x ⊖ u
[chat, vc_chop]    = fl_chop(tempDiff_chop / ww, kk);  % (x ⊖ u) ⊘ w

%% Display results (chopping)
disp('5-digit chopping results:')
fprintf('(a) x ⊕ y  = %g\n', ahat);
fprintf('(b) x ⊖ u  = %g\n', bhat);
fprintf('(c) (x ⊖ u) ⊘ w = %g\n\n', chat);

%% Problem 4: Compute values using 5 digit rounding
% (a) x ⊕ y
[ahat_r, ~] = fl_round(xx + yy, kk);

% (b) x ⊖ u
[bhat_r, ~] = fl_round(xx - uu, kk);

% (c) (x ⊖ u) ⊘ w  (DO round after subtraction, then after division)
[tempDiff_round, ~] = fl_round(xx - uu, kk);    % x ⊖ u (rounded)
[chat_r, vc_round]  = fl_round(tempDiff_round / ww, kk); % (x ⊖ u) ⊘ w

%% Display results (rounding)
disp('5-digit rounding results:')
fprintf('(a) x ⊕ y  = %g\n', ahat_r);
fprintf('(b) x ⊖ u  = %g\n', bhat_r);
fprintf('(c) (x ⊖ u) ⊘ w = %g\n', chat_r);