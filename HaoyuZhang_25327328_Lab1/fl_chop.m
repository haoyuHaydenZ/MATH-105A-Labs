%% Float Chop Funcion
function [fy, vv] = fl_chop(yy, kk)

    if yy == 0
        vv = [zeros(kk,1); 0];
        fy = 0;
        return;
    end

    sgn = sign(yy);
    yy  = abs(yy);

    %Get vector representation (mantissa digits + exponent)
    %Use k-1 digits after the decimal so the mantissa has exactly k digits total
    str_format = sprintf('%%.%de', kk-1);
    mystr_yy   = sprintf(str_format, yy);

    e_indx = strfind(mystr_yy, 'e'); % Get exponent n from the 'e' part (no +1 needed since we used kk-1)
    n      = str2double(mystr_yy(e_indx+1:end));
    vv     = zeros(kk+1,1);
    vv(end)= n;

    % Remove the decimal and read the first k digits of the mantissa
    mantissa = mystr_yy(1:e_indx-1);
    mantissa = strrep(mantissa, '.', '');
    for ii = 1:kk
        vv(ii) = str2double(mantissa(ii));
    end

    pow = kk - 1 - n;
    fy  = sgn * ( floor( yy * 10^pow ) / 10^pow );
end
