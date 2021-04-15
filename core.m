close all; clear all;

[w, fs_rec] = audioread("./recordings/fc_690MHz_fs_10MHz.wav",[10000, 10000*200]);
iq = w(:,1) + j*w(:,2);

% DVB-T parameters
M8K_T = (7/64)*1e-6;    % 8MHz elementary period T = 7/64 microsecs
M8K_FS = 1/M8K_T;       % 8MHz sampling frequency (1/T)
M8K_Su = 8192;          % 8MHz useful part without guard interval (in samples)
M8K_Tu = M8K_Su*M8K_T;  % 8MHz useful part without guard interval (in secs)
M8K_OSg = [2048, 1024, 512, 256]; % 8MHz options for guard interval in samples
M8K_carriers = 6817;    % 8MHz number of carriers in each OFDM symbol [0;6816]
SympF = 68;             % Symbols per frame
FpSF = 4;               % Frames per SuperFrame
load('M8K_ConPilCarr.mat')% Location of continual pilot carriers (8K mode)


% Other params
OPT_RESAMPLE = 1;   % Downsample to M8K_FS
OPT_ADCCorr  = 0;   % Compensate ADC time drift
OPT_BPFilter = 1;   % BP filter only useful subcarriers

BP_margin = 100;    % num of subcarriers additionaly ommited by BP filter to avoid data loss in case of CFO shift (two-sided)
BP_filt_len=500;    % length of BP filter impulse response

bQ = 3;        % Number of Tu memory in blind mode detection
bN = 100e3;    % Number of processed samples during blind mode detection
bNi = 8192;    % Tested Useful Part length

sL = 2;        % Number of observed symbols during finding first symbol
sNi = 8192;    % Useful Part during finding first symbol
sDelta = 0;

aN = M8K_Su*200;    % Number of samples during ADC correction
aWAdd = 200;        % Additional samples added to Gi window length
aGi   = 2048;       % manual Gi set to avoid two-times 1stSymb computation

CFOSymAvg = SympF;      % Number of symbols during CFO integer estimating
CFOMaxIntDev = 200;     % Max multiply of carr spacing to check

DUMMYShift = 0;         % For debugging purposes

%#######################################

figure(1); pwelch(iq(1:2^16),2048,2048-1024,2048,fs_rec,'centered'); title('raw received IQ signal'); 

% DC removal
iq = iq - mean(iq);
figure(2); pwelch(iq(1:2^16),2048,2048-1024,2048,fs_rec,'centered'); title('after DC removal'); 

% Signal gaining
P=2; s = sqrt((length(iq)*P)/(sum(iq.*conj(iq)))); iq=s*iq;

% resample to elementary period T
if OPT_RESAMPLE == 1
    converter = dsp.SampleRateConverter('Bandwidth', ((M8K_carriers+400)/M8K_Su)*M8K_FS, ...
        'InputSampleRate', fs_rec, 'OutputSampleRate', M8K_FS);
    [interp_fact, decim_fact] = getRateChangeFactors(converter);
    iq=converter(iq(1:length(iq)-mod(length(iq),35))); %iq len must be multiple of decimation factor
    %visualizeFilterStages(converter);
end
figure(3); pwelch(iq(1:2^16),2048,2048-1024,2048, M8K_FS,'centered'); title('after resampling to 1/T'); 

% BP filtering - leave used subcarriers only + CFO margin
if OPT_BPFilter == 1
    nyq_freq = floor(M8K_FS/2);
    data_carriers = (((M8K_carriers/M8K_Su)*M8K_FS)/2)/nyq_freq;
    margin = 0.5*(BP_margin * ceil(M8K_FS/M8K_Su))/nyq_freq; % Don't cut-off data carriers when CFO exists
    cutoff_freq=data_carriers+margin;
    bp_h = fir1(BP_filt_len, cutoff_freq, 'bandpass', chebwin(BP_filt_len+1));
    %iq=filtfilt(bp_h,1,iq);
    iq=conv(iq,bp_h); iq=iq(1+floor(length(bp_h)/2):end-floor(length(bp_h)/2));
    figure(4); pwelch(iq(1:2^16),2048,2048-1024,2048, M8K_FS,'centered'); title('after BP filter'); 
end


% ADC correction
if OPT_ADCCorr == 1
    Wsize = aGi+aWAdd;
    [~, adc_acorr] = AutoCorr(iq(1:aN), M8K_Su, bQ*bNi);
    [~, adc_avgacorr] = AvgAutoCorr(adc_acorr, sL, sNi);
    [f_max, f_idx] = max(abs(adc_avgacorr(1:M8K_Su+Wsize)));
    jump_end = floor((length(adc_avgacorr)-(f_idx+Wsize)) / (M8K_Su+aGi));
    peaks=[f_idx];
    for jump=1:jump_end;
        jstart = f_idx+jump*(M8K_Su+aGi)-Wsize/2;
        jend = f_idx+jump*(M8K_Su+aGi)+Wsize/2;
        
        [jf_max, jf_idx] = max(abs(adc_avgacorr(jstart:jend)));
        jf_idx=jf_idx+jstart;
        peaks = [peaks jf_idx];
    end
    
    shouldBe=(length(peaks)-1)*(M8K_Su+aGi);
    reallyIs=peaks(end)-peaks(1);
    step=reallyIs/shouldBe; N=length(iq); new_iq=zeros(1,length([0:step:N-1]));
    new_iq(1,1:length([0:step:N-1]),1)=interp1( [0:N-1], iq, [0:step:N-1], 'spline' );
    iq=new_iq.'; clear new_iq;
    
end

% Blind mode detection
[x_shift, x] = AutoCorr(iq(1:bN), bNi, bQ*bNi);
figure(5); plot(((1:length(x))+x_shift), abs(x)); title('x - AutoCorrelation');

% Highest variation-to-avg ratio = proper (2K/8K) mode
M = ( mean(x.*conj(x)) - mean(abs(x)).*2 ) / mean(abs(x));

% Determine guard interval length
GIpos = islocalmax(abs(x),'MaxNumExtrema', floor(bN/(bNi+max(M8K_OSg))), ...
    'MinSeparation', bNi+min(M8K_OSg)-100); %TODO change bNi to a detected value % slooow 'MinSep'
p=find(GIpos==1); if mod(length(p),2)~=0 p=p(1:end-1); end
p=reshape(p, [2,(length(p)/2)]); p=p(2,:)-p(1,:);
[~,Giidx] = min(abs(M8K_OSg-mean(p))); 

DSGi = M8K_OSg(Giidx);  % Detected Samples per GI
SpSym = M8K_Su + DSGi;  % Samples per symbol (Useful+GI)
SpF = SpSym * SympF;    % Samples per Frame
SpSF = SpF * FpSF;      % Samples per SuperFrame

% Find first sample of the arbitrary symbol
[xf_shift, xf] = AvgAutoCorr(x, sL, sNi);
figure(6); plot((1:length(xf))+xf_shift+x_shift, abs(xf)); title('x̄ - Avg AutoCorrelation');
[~, n1stSymb] = max(xf(1:2*sNi+1)); %first sample of symbol without Gi
n1stSymb=n1stSymb+x_shift+xf_shift+DUMMYShift; %iq index correction


% CFO estimation (3-step approach)
figure(7); plot(abs(10*log10(fftshift(fft(iq(n1stSymb:n1stSymb+M8K_Su)))))); title("CFO"); hold on;
b0 = angle(xf(n1stSymb))/(2*pi); %fraction
idx_CFOiq = n1stSymb:n1stSymb+CFOSymAvg*SpSym-DSGi-1;
iq = iq .* exp(-j*2*pi*(b0/M8K_Tu)*((0:length(iq)-1)/M8K_FS)).';


UsSymbs = getUsefull(iq(idx_CFOiq), M8K_Su, DSGi, 'usefull'); % remove Gi
R=DVBFFT(UsSymbs.', M8K_Su, M8K_carriers).'; %Ceils
p=zeros(1,2*CFOMaxIntDev+1);
for k0=-CFOMaxIntDev:1:CFOMaxIntDev
    kp = M8K_ConPilCarr(9:end-2) + 1; % 200+ left-right possible
    avg_RconjR = mean( R(2,kp+k0) .* conj(R(1,kp+k0)) );
    avg_Rj2 = mean( abs(R(2,kp+k0)).^2 );
    avg_Rj1 = mean( abs(R(1,kp+k0)).^2 );
    p(k0+CFOMaxIntDev+1)=avg_RconjR / sqrt(avg_Rj2 * avg_Rj1);
end
Kidx=find(p==max(p)); K0=Kidx-(CFOMaxIntDev+1); %Integer CFO shift
iq = iq .* exp(-j*2*pi*(K0/M8K_Tu)*((0:length(iq)-1)/M8K_FS)).';% Integer CFO compensation


%step 3 calculate b0 second time, but first retrieve x (Blind mode detection) basing on CFO corrected iq
[s2nd_x_shift, s2nd_x] = AutoCorr(iq(1:bN), M8K_Su, bQ*bNi);
[s2nd_xf_shift, s2nd_xf] = AvgAutoCorr(s2nd_x, sL, M8K_Su);
[~, n1stSymb] = max(s2nd_xf(1:sNi+1)); %first sample of symbol without Gi
n1stSymb=n1stSymb+s2nd_x_shift+s2nd_xf_shift+DUMMYShift; %iq index correction

s2nd_b0 = angle(s2nd_xf(n1stSymb))/(2*pi); %fraction
idx_CFOiq = n1stSymb:n1stSymb+CFOSymAvg*SpSym-DSGi-1;
iq = iq .* exp(-j*2*pi*(s2nd_b0/M8K_Tu)*((0:length(iq)-1)/M8K_FS)).'; % fractional
plot(abs(10*log10(fftshift(fft(iq(n1stSymb:n1stSymb+M8K_Su)))))); legend({'Original', 'After CFO correction'});
figure(8); plot((-CFOMaxIntDev:1:CFOMaxIntDev),abs(p)); title('Integer CFO detection basing on continuous pilots')

UsSymbs = getUsefull(iq(idx_CFOiq), M8K_Su, DSGi, 'usefull'); % remove Gi
sR=DVBFFT(UsSymbs.', M8K_Su, M8K_carriers).'; %Ceils



[idx_testx, testx]=AutoCorr(iq, M8K_Su, bQ*bNi);
[idx_testxf, testxf] = AvgAutoCorr(testx,sL, sNi);
test1stSymb = find(testxf(1:sNi+1)==max(testxf(1:sNi+1))); %first sample of symbol without Gi


