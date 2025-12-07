function [fy, vv] = fl_round(yy, kk)
% This function performs k-digit rounding
% Inputs:
%   yy - real number
%   kk - natural number
% Outputs:
%   fy - f_l_r(yy)
%   vv - [d_1, d_2, ..., d_k, n] finite digit representation of yy

    vv = zeros(kk+1,1);

    str_format = sprintf('%%.%de', kk);
    mystr_yy = sprintf(str_format, yy);

    e_indx = strfind(mystr_yy, 'e');
    nn = str2num(mystr_yy((e_indx+1):end)) + 1;

    yy = yy + (5 * 10^(nn - kk - 1));% Add half ULP for rounding

    [fy, vv] = fl_chop(yy, kk);% Perform chopping after rounding adjustment

end
