%% Q5 Results Tables
A_exact = xx + yy;          % (a) x ⊕ y
B_exact = xx - uu;          % (b) x ⊖ u
C_exact = (xx - uu) / ww;   % (c) (x ⊖ u) ⊘ w

%Chopping estimates already computed above: ahat, bhat, chat
names = {'(a) x ⊕ y'; '(b) x ⊖ u'; '(c) (x ⊖ u) ⊘ w'};

est_chop   = [ahat;  bhat;  chat];
actual_val = [A_exact; B_exact; C_exact];
act_err_c  = est_chop - actual_val;
abs_err_c  = abs(act_err_c);
rel_err_c  = abs_err_c ./ abs(actual_val);

T_chop = table(names, est_chop, actual_val, act_err_c, abs_err_c, rel_err_c, ...
    'VariableNames', {'Operation','Result','ActualValue','ActualError','AbsoluteError','RelativeError'});

disp('--- Results: 5-digit CHOPPING ---');
disp(T_chop);

[tempDiff_round, ~] = fl_round(xx - uu, kk);
[chat_r, ~]         = fl_round(tempDiff_round / ww, kk);

%Round estimates
est_round  = [ahat; bhat; chat_r];
act_err_r  = est_round - actual_val;
abs_err_r  = abs(act_err_r);
rel_err_r  = abs_err_r ./ abs(actual_val);

T_round = table(names, est_round, actual_val, act_err_r, abs_err_r, rel_err_r, ...
    'VariableNames', {'Operation','Result','ActualValue','ActualError','AbsoluteError','RelativeError'});

disp('--- Results: 5-digit ROUNDING ---');
disp(T_round);

writetable(T_chop,'chop_results.csv'); 
writetable(T_round,'round_results.csv');