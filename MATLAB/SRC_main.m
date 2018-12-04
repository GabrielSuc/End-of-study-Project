clc
clear

%This script is meant to provide a Sample Rate Conversion between any
%desired frenquencies. The choice of the filtercan be decided by evaluation
%the differents following cells. The emphasis has been put on FIR
%Parks-Mclellan and IIR Elliptic filters since they both perform best (in 
%terms of MPOS )in their respective category.

%Different specifications for the filter can be chosen such as the ripple  
%in the passband, the attenuation in the stopband and the transition width.
%Moreover, one can also choose how to implement the filters by direct
%implementation, polyphase or multistage.

%Finally, the quality of the conversion can be assessed by using Geoff
%Martin's tool (cosine_sweep_evaluation_tool.m) by "passing" through this
%and his script the filter implementation we desire.

%Gabriel Suc, Bang&Olufsen 2018




%% Frequencies specification

%Type here the desired sampling frequencies that need to be converted

%Input frequency
Fsin = 44.1e3; 
%Output frequency
Fsout = 2*Fsin; %48e3; %

% Resampling Factors L/M
[L,M] = getSRFactors(Fsin,Fsout); % Fsout = (L/M)*Fsin



%%  Lowpass Filter Specifications

% Minimum passband gain i.e. passband ripple 
Rp = 0.1; %in dB

% Maximum stopband gain i.e. stopband attenuation
Rs = 140; %in dB

%Desired Transition Width
%Represent how close we want the passband frequency to be from the stopband
%frequency.
TW = 0.85; %in percent of the stopband freqeuncy;

%Cutoff frequency
%Also defined the stopband frequency i.e. Fstop
Wc = min(pi/M, pi/L);

% Underlying continous-time sampling frequency 
Fe = L*Fsin;

% Equivalent continous-time cutoff frequency
Fc = Wc*(Fe/(2*pi));

%Passpand frequency
Fp = Fc*TW;



%% Sweeps

%Here we generate the cosine sweep (actually a sine now) as our input signal
%This code has be taken from the cosine sweep evaluation tool made by Geoff
%Martin


%We will create several type of signal: for the moment 4 types
%dBFS varies from 0 to -1 and the bit depth from 16 to 24

%cosine_sweeps = zeros(4,1);

fs = Fsin; %_Choose Sampling Frequency
    %_for fs = [44100 48000 88200 96000 176400 192000]
 
    
i = 1; %add counter for filling in signal array
    

for bit_depth = 16 %[16 24], no difference between the two
%        for bitrate = [128 256 320]
    
    audio_channels = 2;         % channels
    f_start = 20;               % Hz
    f_stop = fs/2 * (20000/22050);          % Hz (forces fmax to be 20k @ 44.1 kHz)
    duration_in_seconds = 60;   % seconds
    bitrate = 256;              % AAC bitrate - options: 64, 96, 128, 160, 192, 256, 320 kbps
    for level_dBFS = 0 %-1]    % dB
        
        
        %cosine_sweep{i,1} = struct;
        
        cosine_sweeps(i).fs = fs;
        cosine_sweeps(i).bit_depth = bit_depth;
        cosine_sweeps(i).bitrate = bitrate;
        cosine_sweeps(i).level_dBFS = level_dBFS;

        zero_padding = 0;           % seconds
        dirac_on = 0;               % 1 = on, 0 = off\

        cosine_sweeps(i).sweep = cosine_sweep(audio_channels, f_start, f_stop, duration_in_seconds,... 
                                    level_dBFS, fs, bit_depth, zero_padding, dirac_on);
    
        
       audiowrite(['~/Documents/End-of-study-Project/MATLAB/' num2str(round(cosine_sweeps(i).fs/1000)) 'k_' ...
    num2str(cosine_sweeps(i).bit_depth) '_' num2str(cosine_sweeps(i).level_dBFS) 'dBFS_input.wav'],...
    cosine_sweeps(i).sweep, cosine_sweeps(i).fs, 'BitsperSample', cosine_sweeps(i).bit_depth)                             
         
          

        i = i + 1;
    end
    %end
end

   


       






%% Multistage Implementation

clc

%Take one of the filter above and implement it as a multistage design

%First, determine which combination is the best, for how many stages and
%what kind of combination (in terms of filter: Elliptic, P-M or combination
% of one elliptic + Schuessler filters (minimum-phase filters))

%Regarding the choice of the values of the different stages, it's done by
%hand in getListStages.m. 

[bestPerm,filter_choice,multistage_method] = multi_stage(L,M,Fsin,Fsout,Fp,Rp,Rs); 

Input = input('Choose the input you desire: [1] cosine sweep, [2] impulse, [3] sine: ');

if isempty(Input)
    Input = '1';
end

if Input == 1
    
    %Choose how many samples to plot
    nbr_samples = 10000;
    
    for i = 1:length(cosine_sweeps)
       
        input = cosine_sweeps(i).sweep;

        signal = multistage(Fsin,Fsout,Fp,Rp,Rs,input,bestPerm,filter_choice,multistage_method,nbr_samples,Input);
        signal = signal/max(abs(signal(:))); 


        %Writting the resulting signal as an audio file
        audiowrite(['~/Documents/End-of-study-Project/Sweeps/' num2str(round(Fsout/1000)) 'k_' ...
        num2str(cosine_sweeps(i).bit_depth) '_' num2str(cosine_sweeps(i).level_dBFS) 'dBFS.wav'],...
        signal, floor(cosine_sweeps(i).fs*(L/M)), 'BitsperSample', cosine_sweeps(i).bit_depth)
    
    end
    
    
elseif Input == 2 
    
    %Choose how many samples to plot
    nbr_samples = 10000;
    %Creating Impulse
    t = (-1:(1/Fsin):1)';
    impulse = t==0;
    input = impulse; 
    
    signal = multistage(Fsin,Fsout,Fp,Rp,Rs,input,bestPerm,filter_choice,multistage_method,nbr_samples,Input);
    signal = signal/max(abs(signal(:))); 
    
else
    
    %Choose how many samples to plot
    nbr_samples = 6;
    %Creating Impulse
    t = (0:(1/(Fsin)):1/nbr_samples)';
    f_desired = input('Choose the frequency of the cosine [Hz]: ');
    input = sind(2*pi*f_desired.*t); 
    
    signal = multistage(Fsin,Fsout,Fp,Rp,Rs,input,bestPerm,filter_choice,multistage_method,nbr_samples,Input);
    signal = signal/max(abs(signal(:))); 
    
end 




%% Plots and Saving

figure
%Input and Output

if Input == 1
    plot((0:1/Fsin:(nbr_samples-1)/Fsin),cosine_sweeps(i).sweep(1:nbr_samples))
    hold on
    %Have to compensate for difference in sample rates (to have same
    %ending)
    plot((0:1/Fsout:((nbr_samples-1)/Fsout) + (nbr_samples/Fsin - nbr_samples/Fsout))...
        ,signal(1:(nbr_samples + (nbr_samples/Fsin - nbr_samples/Fsout)*Fsout)))
    legend(['Input at ', num2str(Fsin)], ['Output at ', num2str(Fsout)])
    title('Input and Output signals');
    grid on;
    
elseif Input == 2
    
    plot(t,impulse)
    hold on
    plot(signal)
    legend(['Input at ', num2str(Fsin)], ['Output at ', num2str(Fsout)])
    title('Input and Output signals');
    grid on;
    
else 
    
    plot(t,input)
    hold on
    plot((0:1/(Fsout):1/((Fsout/Fsin)*nbr_samples)), signal(1:length((0:1/(Fsin):1/nbr_samples))))
    legend(['Input at ', num2str(Fsin)], ['Output at ', num2str(Fsout)])
    title('Input and Output signals');
    grid on;
    
end



%snr(signal(:,1), Fsin)







 
