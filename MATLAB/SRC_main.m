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




%%  Frequencies specification

%Type here the desired sampling frequencies that need to be converted

%Input frequency
Fsin = 44.1e3; 
%Output frequency
Fsout =  48e3; %2*Fsin; %

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
TW = 0.95; %in percent of the stopband freqeuncy;

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

%Here we generate the cosine sweep as our input signal
%This code has be taken from the cosine sweep evaluation tool made by Geoff
%Martin


%We will create several type of signal: for the moment 4 types
%dBFS varies from 0 to -1 and the bit depth from 16 to 24

%cosine_sweeps = zeros(4,1);

fs = Fsin; %_Choose Sampling Frequency
    %_for fs = [44100 48000 88200 96000 176400 192000]
 
    
i = 1; %add counter for filling in signal array
    
for bit_depth = [16 24]
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
    
        
%        audiowrite(['~/Documents/End-of-study-Project/Sweeps/' num2str(round(cosine_sweeps(i).fs/1000)) 'k_' ...
%     num2str(cosine_sweeps(i).bit_depth) '_' num2str(cosine_sweeps(i).level_dBFS) 'dBFS_input.wav'],...
%     cosine_sweeps(i).sweep, cosine_sweeps(i).fs, 'BitsperSample', cosine_sweeps(i).bit_depth)                             
                                
        i = i + 1;
    end
    
end

   


       

%% Parks-McClellan Direct & Folded Implementation


%firpmord uses the algorithm suggested in  Rabiner, Lawrence R., and 
%Otto Herrmann. “The Predictability of Certain Optimum Finite-Impulse-Response 
%Digital Filters.” IEEE® Transactions on Circuit Theory. Vol.  20, Number 4, 
%1973, pp. 401–408. This method is inaccurate for band edges close to 
%either 0 or the Nyquist frequency, fs/2.

% cutoff frequencies
f = [Fp Fc];   
% desired amplitudes
a = [1 0];    

% compute deviations
dev = [(10^(Rp/20)-1)/(10^(Rp/20)+1)  10^(-Rs/20)];  

% Find required order to meet the specs
[n_pm_direct,fo,ao,w] = firpmord(f,a,dev,Fe);

% Need to adapt the gain
PM_direct = L*firpm(n_pm_direct,fo,ao,w); %Direct Implementation


%fvtool(PM_direct) %Uncomment to see graphs

%Folded Implentation FIR
PM_folded = dfilt.dfsymfir(PM_direct);

%fvtool(PM_folded)






%% Filtering Part

for i = 1:length(cosine_sweeps)
    
    %Upsampling
    input_uped = upsample(cosine_sweeps(i).sweep,L);
    %Filtering
    input_filtered = filter(PM_folded,input_uped);
    %Downsampling
    signal = downsample(input_filtered,M);
    
    signal = signal/max(signal); %_To prevent data from clipping when writing file
    
    %Writting the resulting signal as an audio file
    audiowrite(['~/Documents/End-of-study-Project/Sweeps/' num2str(round(cosine_sweeps(i).fs/1000)) 'k_' ...
        num2str(cosine_sweeps(i).bit_depth) '_' num2str(cosine_sweeps(i).level_dBFS) 'dBFS.wav'],...
        signal, cosine_sweeps(i).fs, 'BitsperSample', cosine_sweeps(i).bit_depth)
end    
    
%Then we can use Geoff's cosine sweep evaluation tool with tool = 3 and
%process those files

%This part can be not relevant since there's too many coefficients for some
%cases of L and M




%% Elliptic Direct & Cascaded Filters


%firpmord uses the algorithm suggested in Rabiner, Lawrence R., and B. Gold. 
%Theory and Application of Digital Signal Processing. Englewood Cliffs,
%NJ: Prentice-Hall, 1975.


[n_ellip_direct,Wp_ellip_direct] = ellipord(Fp/(Fe/2),Fc/(Fe/2),Rp,Rs);

[z_ellip_direct,p_ellip_direct,k_ellip_direct] = ellip(n_ellip_direct,Rp,Rs,Wp_ellip_direct);

[b_ellip_direct,a_ellip_direct] = ellip(n_ellip_direct,Rp,Rs,Wp_ellip_direct);


Ellip_direct = zp2sos(z_ellip_direct,p_ellip_direct,L*k_ellip_direct); %Have to compensate for he gain

%fvtool(Ellip_direct) %Uncomment to see graphs

% Direct-form II SOS 

Ellip_cascaded = dfilt.df2sos(Ellip_direct);
% 
% Ellip_cascaded.persistentmemory = true;

%fvtool(Ellip_cascaded)

%Other way to implement
% Hlp1 = design(filtSpecs, 'ellip', 'FilterStructure', 'df2sos', ...
%     'MatchExactly', 'both','SOSScaleNorm', 'Linf', ...
%     'SystemObject', true);





%% Filtering Part

for i = 1:length(cosine_sweeps)
    
    %Upsampling
    input_uped = upsample(cosine_sweeps(i).sweep,L);
    %Filtering
    input_filtered = filter(Ellip_cascaded,input_uped); %filter works only with cascaded for some small M and L, sosfilt otherwise
    %Downsampling
    signal = downsample(input_filtered,M);
    
    signal = signal/max(signal); %_To prevent data from clipping when writing file
    
    %Writting the resulting signal as an audio file
    audiowrite(['~/Documents/End-of-study-Project/Sweeps/' num2str(round(cosine_sweeps(i).fs/1000)) 'k_' ...
        num2str(cosine_sweeps(i).bit_depth) '_' num2str(cosine_sweeps(i).level_dBFS) 'dBFS.wav'],...
        signal, cosine_sweeps(i).fs, 'BitsperSample', cosine_sweeps(i).bit_depth)
end      

%Then we can use Geoff's cosine sweep evaluation tool with tool = 3 and
%process those files

%If any problem occurs, go to SNR part down here to have a deeper look into
%what's going on




%% Polystage implementation

%% FIR Polyphase
%I decided to use the same polyphase decomposition given in Russell's
%technique. Even if it is suboptimal it is rather simple to implement and
%still gives a good result in terms of efficiency regarding other
%techniques cf. "Generalized Rational Sampling Rate Conversion Polyphase FIR Filter".


%We need to have the coefficient a and b of the delays

a = 0;
b = 0.1;

while(rem(b,1)~=0)
    a = a + 1;
    b = (L*a - 1)/M;
end   


disp('----------------------- Delay Coefficients ------------------------')
disp( ['a = ', num2str(a), ' and b = ', num2str(b)]);
disp('-------------------------------------------------------------------')

%Now need the LM polyphase components ek of hn

%ATTENTION! POLYPHASE WON'T WORK HERE FOR CERTAIN VALUE OF L AND M.
%THE LENGTH OF THE FILTER HAS TO BE BIGGER THAN L*M ...

ek = myPolyphase(PM_direct,1,L,M,'2');

%Have to filter xin through each branch 

sumBranch = 0;

for k = 1:length(cosine_sweeps)

    for i = (L*M):-1:1 %Starting from the LM-1 branch
        %Creating the other branches before summation
        delayedBy_a = delayseq(cosine_sweeps(k).sweep,(i-1)*a); 
        downsamp = downsample(delayedBy_a,M);
        %Implement polyphase components as folded structures
        filter_polyphase = filter(ek(i,:),1,downsamp);
        upsamp = upsample(filter_polyphase,L); 

        %Sum 

        sumBranch = sumBranch + upsamp;

        if i > 1
            sumBranch = delayseq(sumBranch,-b);
        end
    end
    
    sumBranch = sumBranch/max(sumBranch); %_To prevent data from clipping when writing file
    
    %Writting the resulting signal as an audio file
    audiowrite(['~/Documents/End-of-study-Project/Sweeps/' num2str(round(cosine_sweeps(k).fs/1000)) 'k_' ...
        num2str(cosine_sweeps(k).bit_depth) '_' num2str(cosine_sweeps(k).level_dBFS) 'dBFS.wav'],...
        sumBranch, cosine_sweeps(k).fs, 'BitsperSample', cosine_sweeps(k).bit_depth)
    
end


%% IIR Polyphase

%This part consist of implementing the above filters as a polyphase
%implementation. This can be better explained in "Discrete-Time Signal 
%Processing" by Alan V. Oppenheim& Ronald W. Schafer.

%We split the nonrecursive part (the zeros) and commute it across the
%expander. only if EXPANSION RATE IS HIGHER THAN THE DECIMATION RATE (L > M)



for i = 1:length(cosine_sweeps)
    
    %Polyphase of the nonrecursive part 
    poly_comp = myPolyphase(z_ellip_direct,1,L,M,'1');
    %Filter by the polyphase components
    input_filtered_1st = filter(dfilt.dfsymfir(poly_comp),cosine_sweeps(i).sweep);
    %Upsampling
    input_uped = upsample(input_filtered_1st,L);
    %Filter by the poles as folded structure
    input_filtered_2nd = filter(1,p_ellip_direct,input_uped);
    %Downsampling
    signal = downsample(input_filtered_2nd,M);

    signal = signal/max(signal); %_To prevent data from clipping when writing file
    
    %Writting the resulting signal as an audio file
    audiowrite(['~/Documents/End-of-study-Project/Sweeps/' num2str(round(cosine_sweeps(i).fs/1000)) 'k_' ...
        num2str(cosine_sweeps(i).bit_depth) '_' num2str(cosine_sweeps(i).level_dBFS) 'dBFS.wav'],...
        signal, cosine_sweeps(i).fs, 'BitsperSample', cosine_sweeps(i).bit_depth)
end      

%Then we can use Geoff's cosine sweep evaluation tool with tool = 3 and
%process those files

%If any problem occurs, go to SNR part down here to have a deeper look into
%what's going on




%% Multistage Implementation

clc

%Take one of the filter above and implement it as a multistage design

%First, determine which combination is the best, for how many stages and
%what kind of combination (in terms of filter: Elliptic, P-M or combination
% of one elliptic + Schuessler filters (minimum-phase filters))

%Regarding the choice of the values of the different stages, it's done by
%hand in getListStages.m. 

[bestPerm,manual] = multi_stage(L,M,Fsin,Fsout,Fp,Rp,Rs);

for i = 1:length(cosine_sweeps)

signal = multistage(L,M,Fsin,Fsout,Fp,Rp,Rs,cosine_sweeps(i).sweep,bestPerm,manual);

signal = signal/max(abs(signal(:))); 

%snr(signal(:,1), Fsin)

%Writting the resulting signal as an audio file
audiowrite(['~/Documents/End-of-study-Project/Sweeps/' num2str(round(Fsout/1000)) 'k_' ...
    num2str(cosine_sweeps(i).bit_depth) '_' num2str(cosine_sweeps(i).level_dBFS) 'dBFS.wav'],...
    signal, cosine_sweeps(i).fs, 'BitsperSample', cosine_sweeps(i).bit_depth)


end


