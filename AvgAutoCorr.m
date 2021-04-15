function [x_avg_shift, x_avg] = AvgAutoCorr(x,L,Su)

    xf = zeros(1, length(x)-(L-1)*Su); ns = (L-1)*Su+1;
    for n = ns:length(x)
       xf(n) = sum(x(n-(0:L-1)*Su)) / L;
    end
    xf=xf(ns:end); x_avg_shift=ns-1;
    x_avg=xf;
end