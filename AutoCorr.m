function [x_shift,x] = AutoCorr(iq, Su, WSize)
    % Get iq, shift it back Su samples and sum with oriWSizenal one in WSize
    % samples long window
    N = length(iq);
    r1=iq(1+Su:end); r2=iq;
    n_start = WSize+1; x=zeros(1,N-WSize-Su);
    
    x(n_start) = sum(r1(n_start:-1:n_start-WSize)) .* conj(sum(r2(n_start:-1:n_start-WSize)));
    for n = n_start+1:N-Su
        x(n) = x(n-1) - r1(n-1-WSize).*conj(r2(n-1-WSize)) + r1(n).*conj(r2(n));
    end
    x=x/WSize;
    x=x(n_start:end); x_shift=n_start-1; 
    
end