%% Resampling Caracteristcs

% Sample Rates
Fx = 44.1e3; % Original Frenquency i.e. sample rate of original sequence x[n]
Fy = 48e3; % Sample rate of final sequence y[n]

% Resampling Factors L/M
L = 160; % Fy = (L/M)*Fx
M = 147;


%% Filter Specifications

% Width of transition band
Wt = 5e3; %(Hz)

% Minimum passband gain i.e. passband ripple 
Rp = 1; %(dB)

% Maximum stopband gain i.e. stopband attenuation
Rs = 80;

% Cutoff Frequency
Wc = min(pi/M, pi/L);

% Underlying continous-time sampling frequency 
Fe = L*Fx;

% Equivalent continous-time cutoff frequency
Fc = Wc*(Fe/(2*pi));

% Passband edge
Fp = Fc - Wt/2;

% Stopband edge
Fs = Fc + Wt/2;

% Discrete domain equivalent
Wp = Wc*(Fp/Fc)/pi;
Ws = Wc*(Fs/Fc)/pi;



%% Typical Sinc Filter

N = 24*L;
%h = fir1(1608,1/M,kaiser(1609,7.8573)); %cf. Kaiser documentation with side lobes attenuation of 80 dB
%h = L*h;

%H = dfilt.df1sos(h);
%Lacasc = cascade(h);
%fvtool(h)


%% Elliptical Filter


[n,Wp2] = ellipord(Wp,Ws,Rp,Rs);

[b1,a1] = ellip(n,Rp,Rs,Wp); %Second parameter = ripple in the passband (dB) -> cf. threshold articles 
                            %Play a role in the number of samples involved in the Group delay (Amplitude) 

% Second order sections filter
sos = tf2sos(b1,a1); % L-by-6 Matrix, 2*L gives order of the filter
%fvtool(sos);



% ------------- Other Design Method ------------------


filtSpecs = fdesign.lowpass('Fp,Fst,Ap,Ast',Wp,Ws,Rp,Rs);
%filtSpecs = fdesign.lowpass('N,F3db,Ap,Ast',n,Wc,Rp,Rs);


% -------------- Structure Camparison -------------------

% Direct-form I SOS 

Hlp1 = design(filtSpecs, 'ellip', 'FilterStructure', 'df1sos', ...
    'MatchExactly', 'both','SOSScaleNorm', 'Linf', ...
    'SystemObject', true);

fvtool(Hlp1); 

cost(Hlp1)
info(Hlp1)


% Direct-form II SOS  (is polyphase)

Hlp2 = design(filtSpecs, 'ellip', 'FilterStructure', 'df2sos', ...
    'MatchExactly', 'both','SOSScaleNorm', 'Linf', ...
    'SystemObject', true);
  
fvtool(Hlp2);   

cost(Hlp2) 
info(Hlp2)

% Direct-form I transposed SOS df1tsos 

Hlp3 = design(filtSpecs, 'ellip', 'FilterStructure', 'df1tsos', ...
    'MatchExactly', 'both','SOSScaleNorm', 'Linf', ...
    'SystemObject', true);

fvtool(Hlp3);   

cost(Hlp3)
info(Hlp3)

% Direct-form II transposed SOS 


Hlp4 = design(filtSpecs, 'ellip', 'FilterStructure', 'df2tsos', ...
    'MatchExactly', 'both','SOSScaleNorm', 'Linf', ...
    'SystemObject', true);

fvtool(Hlp4);

cost(Hlp4)
info(Hlp4)

% Cascade minimum-multiplier allpass 

Hlp5 = design(filtSpecs, 'ellip', 'FilterStructure', 'cascadeallpass', ...
    'MatchExactly', 'both','SOSScaleNorm', 'Linf', ...
    'SystemObject', true);

fvtool(Hlp5);  

cost(Hlp5)
info(Hlp5)

% Cascade wave digital filter allpass  

Hlp6 = design(filtSpecs, 'ellip', 'FilterStructure','cascadewdfallpass', ...
    'MatchExactly', 'both','SOSScaleNorm', 'Linf', ...
    'SystemObject', true);

fvtool(Hlp6);   

cost(Hlp6) 
info(Hlp6)



% Impulse response
[hiir,t_samp_iir] = impz(Hlp1);
hiir = L*hiir;  % Passband gain is L






%% Parks-McClellan Filter


%firpmord uses the algorithm suggested in [1]. This method is inaccurate 
%for band edges close to either 0 or the Nyquist frequency, fs/2.

f = [(Fp*Wc)/pi (Fsc*Wc)/pi];   % cutoff frequencies
a2 = [1 0];         % desired amplitudes, linear -3 dB 
% compute deviations
dev = [(10^(Rp/20)-1)/(10^(Rp/20)+1)  10^(-Rs/20)];  

% Find required order to meet specs
[n2,fo,ao,w] = firpmord(f,a2,dev,Fx);

%Freq response
b2 = firpm(n2,fo,ao,w);

% Direct-form FIR
Hd1 = dfilt.dffir(b2);

fvtool(Hd1);

cost(Hd1)
info(Hd1)

% Direct-form FIR transposed
Hd2 = dfilt.dffirt(b2);

fvtool(Hd2);

cost(Hd2)
info(Hd2)

% Direct-form symmetric FIR 
Hd3 = dfilt.dfsymfir(b2);

fvtool(Hd3);

cost(Hd3)
info(Hd3)

% Direct-form antisymmetric FIR 
Hd4 = dfilt.dfasymfir(b2);

fvtool(Hd4);

cost(Hd4)
info(Hd4)

% Overlap-add FIR

% The resulting number of FFT points = (filter length + the block length - 1). 
%The filter is most efficient if the number of FFT points is a power of 2. 

Hd5 = dfilt.fftfir(b2,24); %block length 

fvtool(Hd5);

% cost(Hd5)     
% info(Hd5)        Not supported



% Impulse response
[hfir,t_samp] = impz(Hd1);
hfir = L*hfir;  % Passband gain is L

% ----------- Using DSP Toolbox ----------------

% h = fdesign.lowpass('n,fp,fst', n2, Wp, Ws);
% 
% Hd = design(h, 'equiripple', ...
%     'StopbandShape', 'flat', ...
%     'SystemObject', true);
% 
% fvtool(Hd);

% Can't control the value of Rs principally



%% Multistage Proakis

%Create order of desired filter according to example of multistage decomposition 
%of chapter 11 of Proakis'book


%Defining the multistages (Decimator)
D1 = 21; D2 = 7;

%Filters' parameters
%Filter 1
F1 = Fx/D1;
Fsc1 = Fs - F1;

deltaf1 = (Fsc1 -  Fp)/Fx;

%Filter 2
F2 = F1/D2;
Fsc2 = Fs - F2;

deltaf2 = (Fsc2 -  Fp)/F1;

%Ripples
delta11 = (10^(-Rp))/2;
delta21 = (10^(-Rs));

%Orders
ord1 = round((-10*log10(delta11*delta21)-13)/(14.6*deltaf1)+1);
ord2 = round((-10*log10(delta11*delta21)-13)/(14.6*deltaf2)+1);


%Designing the filters
f1 = [0; 2*(F1*Wc)/(pi*Fx); 2*(Fsc1*Wc)/(pi*Fx); 1];
f2 = [0; 2*(F2*Wc)/(pi*F1); 2*(Fsc2*Wc)/(pi*F1); 1];

b01 = firpm(ord1,f1,ao);
b02 = firpm(ord2,f2,ao);

H01 = dfilt.fftfir(b01);
H02 = dfilt.fftfir(b02);

%Cascading them

H12 = cascade(H01,H02);

%Graphs
fvtool(H01)
fvtool(H02)
fvtool(H12)

%Cost
cost(H12)

% Impulse response
[hproakis,t_proakis] = impz(H12);
hproakis = L*hproakis;  % Passband gain is L
%% Multistage MIT

% 4 stages

L1 = 2; L2 = 2; L3 = 2; L4 = 20;

%Filter 1
Fp1 = Fp/L1; Fs1 = Fs/L1; %(cf. thesis MIT)

f1 = [(Fp1*Wc)/pi (Fs1*Wc)/pi];   % cutoff frequencies
a01 = [1 0];         % desired amplitudes, linear -3 dB 
% compute deviations
dev1 = [(10^(Rp/20)-1)/(10^(Rp/20)+1), 1];

% Find required order to meet specs
[n01,fo01,ao01,w01] = firpmord(f1,a01,dev1,Fx);

%Freq response
b01 = firpm(n01,fo01,ao01,w01);

H1 = dfilt.fftfir(b01);

%Filter 2
Fp2 = Fp/(L1*L2); Fs2 = (Fx/(L1*L2))*(L1-(Fs/Fx));

f2 = [(Fp2*Wc)/pi (Fs2*Wc)/pi];   % cutoff frequencies
a02 = [1 0];         % desired amplitudes, linear -3 dB 

% Find required order to meet specs
[n02,fo02,ao02,w02] = firpmord(f2,a02,dev1,Fx/L1);

%Need to adapt frequency
F1 = Fx/L1;

%Freq response
b02 = firpm(n02,fo02,ao02,w02);

H2 = dfilt.fftfir(b02);

%Filter 3
Fp3 = Fp/(L1*L2*L3); Fs3 = (Fx/(L1*L2*L3))*(L1*L2-(Fs/Fx));

f3 = [(Fp3*Wc)/pi (Fs3*Wc)/pi];   % cutoff frequencies
a03 = [1 0];         % desired amplitudes, linear -3 dB 


% Find required order to meet specs
[n03,fo03,ao03,w03] = firpmord(f3,a03,dev1,F1/L2);

%Need to adapt frequency
F2 = F1/L2;

%Freq response
b03 = firpm(n03,fo03,ao03,w03);

H3 = dfilt.fftfir(b03);

%Filter 4
Fp4 = Fp/(L); Fs4 = (Fx/L)*(L1*L2*L3-(Fs/Fx));

f4 = [(Fp4*Wc)/pi (Fs4*Wc)/pi];   % cutoff frequencies
a04 = [1 0];         % desired amplitudes, linear -3 dB 


% Find required order to meet specs
[n04,fo04,ao04,w04] = firpmord(f4,a04,dev1,F2/L3);

%Freq response
b04 = firpm(n04,fo04,ao04,w04);

H4 = dfilt.fftfir(b04);
 
%Cascade
Hcasc = cascade(H1,H2,H3,H4); 

fvtool(Hcasc)

cost(Hcasc)
info(Hcasc)

[hcasc,t_casc] = impz(Hcasc);
hcasc = L*hcasc;  % Passband gain is L

%% Overlap-save Method

%  nfft in the FFT size parameter as a power of two value greater (typically much greater) than n+1
q = length(b2)/2;
nfft = 2^(q+1);
lenH = length(b2);

nr=(length(yc(:,1)))/nfft;

for i = 1:nr
Hd6(i,:) = ifft(fft(yc(i:i+lenH +(-1),1),nfft) .* fft(b2,nfft));
end

%Hd6 = Overlap_Save_Method(yc(:,1)',b2,4096); % One Channel

%% SRC (Using resampling.m)

% Loading a music file
[y,Fsamp] = audioread('Katty_Perry.mp3'); %Katty_Perry.mp3, Muddy_waters.aiff

x = 30;

yc = y(1:x*Fsamp,:); %take x times Fs sec of original signal 

% Resampling with a rational factor
t = 1:x*Fy;

%rs = resample(yc(5062:6062,1),t(5062:6062), Fs, 160,147,'pchip');

rs = resample(yc,L,M);

% Plot comparison
% PLOT FOR ONLY ONE CHANNEL (LEFT)
subplot(3,1,1)  
plot(t(1:x*Fsamp), yc(:,1),'*')
title('original signal')

subplot(3,1,2)
plot(t(1:x*Fsamp), yc(:,1),'*',t,rs(:,1),'o')
xlabel('Samples')
ylabel('Signal')
legend('Original','Resampled','Location','NorthWest')
title('Resampling with resample.m')

% ------------- METHOD 2 ---------------

% up = upsample(yc,L);
% 
% yc_filtered = conv2(up,sos,'same');
% 
% down = downsample(yc_filtered,M);

n_samp= 0:44099;  % 10240 samples, 0.0232 seconds long

y_res = upfirdn(y,hproakis,L,M);

subplot(3,1,3)
plot(n_samp/Fs,y(1:44100),'*'); 
hold on 
plot(n_samp/(Fs*L/M),y_res(1:44100),'o'); 
hold off;
xlabel('Time (sec)');
ylabel('Signal value');
legend('Original','Resampled')
title('Resampling with upfirdn')


sound(y_res(1:x*Fs),Fy)



%% SRC Using DSP Toolbox

%Using Previous Park-McClellan filter

SRC= dsp.FIRRateConverter(L,M,hiir); %Test with hfir, hiir

cost(SRC)

% SRC_polyphase = polyphase(SRC);
% 
% polyphase(SRC)

input = yc;
output = SRC(input);

plot(n_samp/Fs,y(1:44100),'*'); 
hold on 
plot(n_samp/(Fs*L/M),output(1:44100),'o'); 
hold off;
xlabel('Time (sec)');
ylabel('Signal value');
legend('Original','Resampled')
title('Resampling with ds.FIRRC')

%% Farrow SRC

farrowSampRateConv_3rd = dsp.FarrowRateConverter('InputSampleRate',Fx, ...
    'OutputSampleRate',Fy,'PolynomialOrder',3);

cost(farrowSampRateConv_3rd)

input = yc;
output = farrowSampRateConv_3rd(input);


plot(n_samp/Fs,y(1:44100),'*'); 
hold on 
plot(n_samp/(Fs*L/M),output(1:44100),'o'); 
hold off;
xlabel('Time (sec)');
ylabel('Signal value');
legend('Original','Resampled')
title('Resampling with ds.FarrowRC')

%sound(output(1:30*Fs),Fy)

%% Classic SRC (Multistage, Polyphase)



SRC_classic = dsp.SampleRateConverter('Bandwidth',22.05e3,'InputSampleRate',Fx,'OutputSampleRate',Fy);

visualizeFilterStages(SRC_classic)

cost(SRC_classic)

input = yc;
output = SRC_classic(input);


plot(n_samp/Fs,y(1:44100),'*'); 
hold on 
plot(n_samp/(Fs*L/M),output(1:44100),'o'); 
hold off;
xlabel('Time (sec)');
ylabel('Signal value');
legend('Original','Resampled')
title('Resampling with ds.SRC')

%%
filters_classic_SRC = getFilters(SRC_classic);

% Impulse response
[hclassicSRC,t_samp] = impz(filters_classic_SRC);

%--------------- Tryna recreating the filter used -------------

%THIS IS IT
fsrc = fdesign.rsrc(L,M,'Nyquist',L,Fc,80,44100*160);
hsrc = kaiserwin(fsrc);

[btest,atest] = tf(SRC);
 
testmypoly = myPolyphase(btest,atest,L,M,'1');
comparepoly = polyphase(SRC);

% fvtool(hclassicSRC)
% fvtool(hsrc)
% 
cost(hsrc)
info(hsrc)

cost(filters_classic_SRC)
info(filters_classic_SRC)