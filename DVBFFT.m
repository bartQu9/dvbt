% Get NxM where N - length of each IQ vect ; M - IQ vectors
function ret = DVBFFT(iq, Nfft, Ncarr)
    [N,M] = size(iq);
    if N ~= Nfft error('length(iq) != Nfft'); end
    
    X = fft(iq);
    X_sh = zeros(N,M);
    for m=1:M
        X_sh(:,m) = fftshift(X(:,m));
    end
    S = X_sh(Nfft/2-floor(Ncarr/2)+1 : Nfft/2+floor(Ncarr/2)+1,:);
    ret = S;
end