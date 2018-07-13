clear all; 
close all;
clc;
%% Resampling Caracteristcs

% Sample Rates
Fsin = 44.1e3; % Original Frenquency i.e. sample rate of original sequence x[n]
Fsout = 48e3; % Sample rate of final sequence y[n]

% Resampling Factors L/M
[L,M] = getSRFactors(Fsin,Fsout); % Fy = (L/M)*Fx

%% Filter Specifications


%Choose percentage of bandwidth desired
BW = 0.80;

% Width of transition band
Wt = (2*1102.5); %(Hz)

% Minimum passband gain i.e. passband ripple 
Rp = 0.1; %(dB)

% Maximum stopband gain i.e. stopband attenuation
Rs = 50; %dB

% Cutoff Frequency
Wc = min(pi/M, pi/L);
%Same frequency but not normalized
Fcut = min((Fsin/2)/L,(Fsin/2)/M);

% Underlying continous-time sampling frequency 
Fe = L*Fsin;

% Equivalent continous-time cutoff frequency
Fc = Wc*(Fe/(2*pi));

% Passband edge
Fp = Fc*BW;%Fc - Wt/2;
%Same frequency but not normalized
Fpass = Fcut*BW;

% Stopband edge
Fsc = Fc + Wt/2;

% Discrete domain equivalent
Wp = Wc*(Fp/Fc);%/pi;
Wsc = Wc*(Fsc/Fc);%/pi;

%Specifications for designing filters
filtSpecs = fdesign.lowpass('Fp,Fst,Ap,Ast',Wp,Wsc,Rp,Rs);

%% MULTISTAGE
%Multi-stage

clc;

%[Order, Rt] = multi_stage(L,M,Fsin,Fsout,Fp,Rp,Rs);


%Order_turek = multistage_turek(L,M,Fx,Fp,Fsc,Rp,Rs);


%% FILTERS
%% Raised Cosine with Kaiser Window (Matlab's SRC)

global h_sinc

% fsrc = fdesign.rsrc(L,M,'Nyquist',L,Wc-Wp,Rs,44100); %TW = Wc-Wp
% h_raisedco = kaiserwin(fsrc);%design(fsrc,'SystemObject',true);


%Too long to compute, let's try sinc filter with kaiser window


f3 = [Wp/pi Wc/pi];   % cutoff frequencies
a3 = [1 0];         % desired amplitudes

% compute deviations
dev = [(10^(Rp/20)-1)/(10^(Rp/20)+1)  10^(-Rs/20)]; 

%Let's take the same order we'll have for a Parks-McClellan design
[n3,fo3,ao3,wo3] = firpmord(f3,a3,dev);

h_sinc = firls(n3,fo3,ao3,wo3); %kaiser(n3 + 1,7.8573)); 
h_sinc = L*h_sinc;





fvtool(h_sinc)

%% Elliptic Filters

global Hlp1

% Direct-form II SOS 

% Hlp1 = design(filtSpecs, 'ellip', 'FilterStructure', 'df1sos', ...
%     'MatchExactly', 'both','SOSScaleNorm', 'Linf', ...
%     'SystemObject', true);

[n1,Wp1] = ellipord(Fpass/(Fsin/2),Fcut/(Fsin/2),Rp,Rs);

[z1,p1,k1] = ellip(n1,Rp,Rs,Wp1);

[b1,a1] = ellip(n1,Rp,Rs,Wp1);


Hlp1 = zp2sos(z1,p1,L*k1); %Have to compensate for he gain
%get rid of dfilt.df2sos if we want to use fft, snr and Geoff's analysis

fvtool(Hlp1)

%%
clc;
%Russell

Np = length(p1); %Number of poles

%Have to choose Nl and Nm -> put a munual stuff directly inside function 
%The gain in efficiency is the same regarless of how many poles of H(z)
%areassigned to Hl(z) ans how many bto Hm(z). However, complex-conjugate
%pole pairs should not be separated since this would make the filtr
%coefficients complex. Therefore, choose Nl even and Nm odd if length
%filter odd.

if rem(Np,2) ~=0
    
Nl = input('Choose even number of pole Nl to be assigned to coefficient L: ');

    if isempty(Nl)
          error('Choose Nl');
    
    elseif Nl > Np
        error('Nl has to be smaller than Np')
    
    elseif rem(Nl,2) ~= 0
        error('Nl has to be even')
    end

else
   
    if isempty(Nl)
          error('Choose Nl');
    
    elseif Nl > Np
        error('Nl has to be smaller than Np')
    
    elseif rem(Nl,2) ~= 0
        error('Nl has to be even')
    end
    
end

Nm = Np - Nl;
    
    
xin = [1,1,1,1,1,1,1,1,1,1];

%output_russel = russell(z1,p1,k1,L,M,Nl,Nm,xin);
 

%% Parks-McClellan

global b2

%firpmord uses the algorithm suggested in [1]. This method is inaccurate 
%for band edges close to either 0 or the Nyquist frequency, fs/2.

f = [Fpass Fcut];   % cutoff frequencies
a2 = [1 0];         % desired amplitudes, linear -3 dB 
% compute deviations
dev = [(10^(Rp/20)-1)/(10^(Rp/20)+1)  10^(-Rs/20)];  

% Find required order to meet specs
[n2,fo,ao,w] = firpmord(f,a2,dev,Fsin);

%Freq response
% Need to adapt the gain
b2 = L*firpm(n2,fo,ao,w);


%fvtool(b2)

%global Hd1;

% Direct-form FIR
%Hd1 = dfilt.dffir(b2);



%fvtool(Hd1)









%% Filtering

%FFT Sine

F = 18000;                %Frequency of the sine
Fs = Fsin;                % Sampling frequency                    
T = 1/Fs;               % Sampling period       
len = 44100;               % Length of signal. Has to be bigger than Fs
t = (0:len-1)*T;          % Time vector
 

input_src = 1*sin(2*pi*F*t);%*(0:n-1)); % Original signal %/Fx * n





FFT = fft(input_src);
Y = abs(FFT/len); %length(input_src)

P1 = Y(1:len/2+1);
P1(2:end-1) = 2*P1(2:end-1); %Compensate amplitude division

f = Fsin*(0:len/2)/len;  
plot(f,P1);

%% FFT Upsampled Sine


input_uped = upsample(input_src,L);

FFT = fft(input_uped);
Y = abs(FFT/length(input_uped)); %length(input_src)

P1 = Y(1:length(input_uped)/2+1);
P1(2:end-1) = L*2*P1(2:end-1); %Compensate amplitude division

f = Fsin*L*(0:length(input_uped)/2)/length(input_uped);  

plot(f,P1);


%% FFT Filtered Upsampled Sine 

%Filtering by FIR
%input_filtered = filter(b2,1,input_uped);

%Filtering by IIR
input_filtered = sosfilt(Hlp1,input_uped);


FFT = fft(input_filtered);
Y = abs(FFT/length(input_filtered)); %length(input_src)

P1 = Y(1:length(input_filtered)/2+1);
P1(2:end-1) = 2*P1(2:end-1); %Compensate amplitude division

f = Fsin*L*(0:length(input_filtered)/2)/length(input_filtered);  
plot(f,P1);


%% FFT Downsampled Filtered Upsampled Sine


output_src = downsample(input_filtered,M);


FFT = fft(output_src);
Y = abs(FFT/length(output_src)); %length(input_src)

P1 = Y(1:length(output_src)/2+1);
P1(2:end-1) = 2*P1(2:end-1); %Compensate amplitude division

f = Fsout*(0:length(output_src)/2)/length(output_src); 

plot(f,P1);


%% SNR

% snr(output_src, Fsout);
% 
% 
% newinput = [input_src, zeros(1,length(output_src)-length(input_src))];
% 
% xx=mean(abs((fft(newinput))).^2);
% yy=mean((abs(fft(newinput-dist))).^2);
% mysnr=10*log10(xx/yy);

%%

snr(output_src, Fsin);