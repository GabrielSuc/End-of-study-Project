close all
clear all

% SRC_Comparison

% Tools used to provide a basic introductory analysis of an audio "black
% box". The script can be used to create a cosine sweep that is saved as a
% .wav file and played through a DUT. The output of the DUT is captured and
% analysed with the script, producing a graphic result for basic analysis.
%
% Note that this is a VERY basic test and should not be used to make a go/no-go
% decision about the usability of a device. It can only provide a "no-go" -
% or indications of simple errors in an audio transmission system.
%
% geoff martin, Bang & Olufsen 2017


% tool selection
%
%   1 = generate a cosine wave sweep for basic test
%   2 = generate a file with a click for latency testing
%   3 = analyse the output of the cosine wave sweep test
%   4 = analyse the results of the latency test


tool = 3;



%Everything marked as %_ ... is a modification of the original file made by Geoff 



%%


%%% Cosine sweep generator %%%

if tool == 1
    
    file_type = 1;   % 1=wav, 2=flac, 3 = AAC
    
    fs = 44100; %_Choose Sampling Frequency
    %_for fs = [44100 48000 88200 96000 176400 192000]
        for bit_depth = [16 24]
%        for bitrate = [128 256 320]
            disp('Cosine wave generator')
            
            audio_channels = 2;         % channels
            f_start = 20;               % Hz
            f_stop = fs/2 * (20000/22050);          % Hz (forces fmax to be 20k @ 44.1 kHz)
            duration_in_seconds = 60;   % seconds
            bitrate = 256;              % AAC bitrate - options: 64, 96, 128, 160, 192, 256, 320 kbps
            for level_dBFS = [0 -1]    % dB
                
                fs
                bit_depth
                bitrate
                level_dBFS
                
                zero_padding = 0;           % seconds
                dirac_on = 0;               % 1 = on, 0 = off\
                
                
                h = getFilter('2'); %_Get Filter from SRC_Comparison
                
                %Sample rate conversion
                
                %_
                input_uped = upsample(cosine_sweep(audio_channels, f_start, f_stop, duration_in_seconds,... 
                level_dBFS, fs, bit_depth, zero_padding, dirac_on),L);

                input_filtered = sosfilt(h,input_uped);

                signal = downsample(input_filtered,M);
      

                signal = signal/max(signal); %_To prevent data from clipping when writing file
                %_
                
               
                
                switch file_type
                    case 1
                    audiowrite([num2str(round(fs/1000)) 'k_' num2str(bit_depth) '_' num2str(level_dBFS) 'dBFS.wav'], signal, fs, 'BitsperSample', bit_depth)
                    case 2                    
                    audiowrite([num2str(round(fs/1000)) 'k_' num2str(bit_depth) '_' num2str(level_dBFS) 'dBFS.flac'], signal, fs, 'BitsperSample', bit_depth)
                    case 3                   
                    audiowrite([num2str(round(fs/1000)) 'k_' num2str(bitrate) 'kbps_' num2str(level_dBFS) 'dBFS.m4a'], signal, fs, 'BitRate', bitrate)
                end
                
            end
            
        end
    %_end
    
    
    
    
    
%%% Click generator %%%%
    
    
elseif tool == 2
    
    disp('Click generator for latency testing')
    
    audio_channels = 2;         % channels
    duration_in_seconds = 30;   % seconds
    level_dBFS = -20;           % dB
    fs = 44100;                 % Hz
    bit_depth = 16;             % bits
    
    signal = click_generator(audio_channels, duration_in_seconds, level_dBFS, fs, bit_depth);
    
    audiowrite(['click_track' num2str(duration_in_seconds) 'sec_' num2str(level_dBFS) 'dBFS_' num2str(fs) '_' num2str(bit_depth) '.wav'], signal, fs, 'BitsperSample', bit_depth)
    
        
    

%%% Cosine Sweep Analysis Section %%%%
    
elseif tool == 3

    % sidebar: to detect wideband noise (e.g. tpdf dither) thesne should 
    % be set to approximately 32 dB below the noise floor value
    
    
    % one file, or a directory of files?
    manual = input('1 to analyse a single file, 0 to analyse all audio files in a directory [1]:');
    if isempty(manual)
          manual = 1;
    end
       

    disp('Cosine wave analyser')

    % dynamic range of the analyis in dB.
    % Note that this is the range downwards
    % from the maximum signal - not dB FS!
    analysis_snr = input('Dynamic range of analysis in dB [100]:');
    if isempty(analysis_snr)
          analysis_snr = 100;
    end

    stepsize = 5;                   % difference in dB between contour lines

    
    
    upper_frequency_limit = input('Maximum analysis frequency (set to 0 for Nyquist) [0]:');
    if isempty(upper_frequency_limit)
          upper_frequency_limit = 0;
    end
    
    sweep_length_in_seconds = 0;   % length of audio signal to analyse. Set to 0 for automatic
    trigger_threshold_in_dB = -20;  % level below peak value that signals the start of the sweep
    
    
    
    
    
    
    if manual == 1
        [filename.name, pathname] = uigetfile({'*'}); % ask for and load the file to analyse
        
        filelist(1).name = filename.name;
        
    else
        
        pathname = uigetdir; % ask for and load the file to analyse
        
        cd(pathname)
        
        filelist = dir;     % get the list of filenames
        
        filelist(1:2) = []; % delete the first three hidden files for Mac/ two for linux
        
        m=1;
        
        for n=1:length(filelist)
            if filelist(n).name(end-3:end) == '.wav'
                filename(m).name = filelist(n).name;
                m=m+1;
            elseif filelist(n).name(end-4:end) == '.flac'
                filename(m).name = filelist(n).name;
                m=m+1;
            elseif filelist(n).name(end-3:end) == '.m4a'
                filename(m).name = filelist(n).name;
                m=m+1;
           end
        end
        
    end
    
    % go through the files, one by one
    
    for n=1:length(filelist)
        
        filename_ref = filename(n).name;
        
        
        
        
        
        
        if filename_ref(end-3:end) == '.raw'
            fs = 88200;
            nchannel = 2;
            nbits = 32;
            y = rawaudioread(fullfile(pathname, filename_ref), nchannel, nbits);
        elseif filename_ref(end-3:end) == '.pcm'
            fs = 176400;
            nchannel = 2;
            nbits = 32;
            y = rawaudioread(fullfile(pathname, filename_ref), nchannel, nbits);
        elseif filename_ref(end-3:end) == '.wav' | '.m4a'
            [y, fs] = audioread(fullfile(pathname, filename_ref));
            fs;
        elseif filename_ref(end-4:end) == '.flac'
            [y, fs] = audioread(fullfile(pathname, filename_ref));
            fs;
        end
        
        
        
        
        
        % check if upper_frequency_limit == 0, set to fs/2 if true
        if upper_frequency_limit == 0
            upper_frequency_limit = fs/2;
        end
        
        
        
        
        % clean up the silence at either end of the capture
        y_peak = 20*log10(max(abs(y)));                     % find the peak value in dB FS (sort of....)
        y_start_threshold_linear = 10.^((y_peak + trigger_threshold_in_dB) / 20);
        
        start_sample = min(find(abs(y) > y_start_threshold_linear));
        
        
        
        
        if sweep_length_in_seconds == 0
            y = y(start_sample : end, :);
        else
            if start_sample + sweep_length_in_seconds*fs > length(y)
                y = y(start_sample : end, :);
            else
                y = y(start_sample : start_sample + sweep_length_in_seconds*fs, :);
            end
        end
        
        
        
        
        
        for snr_n = 1:length(analysis_snr) % loop for different SNR values in the analysis
            
            cos_sweep_analyse_and_plot(y, analysis_snr(snr_n), stepsize, fs, upper_frequency_limit, pathname, filename_ref)
            
        end
        
    end
    
 
    close all
    
    
    
    
    
    
%%% Latency Test Analysis Section %%%%
    
elseif tool == 4
    
    savefile_on = 1;
    
    
    detection_threshold_in_dB_rel_peak = -10;
    
    
    [filename, pathname] = uigetfile({'*'}) % ask for and load the file to analyse
    
    filename
    
    
    if filename(end-3:end) == '.raw'
        fs = 44100;
        nchannel = 2
        nbits = 32
        y = rawaudioread(fullfile(pathname, filename), nchannel, nbits);
    elseif filename(end-3:end) == '.wav'
        [y, fs] = audioread(fullfile(pathname, filename));
        fs
    end
    
    
    filename_pdf = filename;    % match the filename of the .wav input and the .pdf output
    filename_pdf(end-3:end) = [];   % strip the suffix from the audio file's name
    filename_pdf = [filename_pdf '_latency.pdf']     % add info about snr to the pdf filename
    
    
    
    
    % check that the file has 2 channels to compare with each other
    
    
    
    
    
    
    
    % find the peak values of the two channels
    
    
    detection_threshold_lin = 10^(detection_threshold_in_dB_rel_peak/20);
    
    y = y.^2;   % convert to "power" and make abs
    
    peaks = max(y);
    
    clicks_1 = find(y(:,1) > peaks(:,1) .* detection_threshold_lin);
    clicks_2 = find(y(:,2) > peaks(:,2) .* detection_threshold_lin);
    
    clicks_1 = clicks_1(find(diff(clicks_1) ~= 1));
    clicks_2 = clicks_2(find(diff(clicks_2) ~= 1));
    
    
    % find the shortest of the two vectors
    if length(clicks_1) < length(clicks_2)
        clicks_2(length(clicks_1)+1:end) = [];
    else
        clicks_1(length(clicks_2)+1:end) = [];
    end
    
    
    t = clicks_2 / fs;
    
    figure
    plot(t, (clicks_2 - clicks_1) ./ fs * 1000000, 'k');    % plots difference in microseconds
    title(filename_pdf(1:end-4), 'Interpreter', 'none');
    ylabel('Time difference (Microseconds)')
    xlabel('Elapsed Time (Seconds)')
    grid on
    temp = axis;
    line([0 temp(2)], [100 100])
    
    
    if savefile_on == 1
        set(gcf, 'PaperPosition', [0 0 21 14]);    %Position plot at left hand corner.
        set(gcf, 'PaperSize', [20 14]);              %Set the paper width and height.
        saveas(gcf, fullfile(pathname, filename_pdf), 'pdf') %Save figure
    end
    
end



%%




function cos_sweep_analyse_and_plot(y, snr, stepsize, fs, upper_frequency_limit, pathname, filename)

channels = 1;    % audio channel number(s). Set to [1 2] for stereo


savefile_on = 1;    % set to 0 to just plot the result without saving a pdf




windowlength = 2^(14+round(fs/44100))    %balancing frequency resolution (high number) vs time resolution (low number)
overlap = windowlength / 4;
snr_step = stepsize;


% check if upper_frequency_limit == 0, set to fs/2 if true
if upper_frequency_limit == 0
    upper_frequency_limit = fs/2;
end


for n=channels   % loop  for the number of audio channels
    [s, f, t] = spectrogram(y(:,n),blackmanharris(windowlength/2), overlap, windowlength, fs, 'yaxis', 'MinThreshold', -1*snr);
    s_dB = 20*log10(abs(s));
    
    
    %    f_in = logspace(log10(10), log10(upper_frequency_limit), length(t));
    
    
    max(max(s_dB));
    min(min(s_dB));
    
    s_dB = s_dB - max(max(s_dB));       % make the analysis relative to the peak in the signal
    
    
    filename_pdf = filename;    % match the filename of the .wav input and the .pdf output
    filename_pdf(end-3:end) = [];   % strip the suffix from the audio file's name
    filename_pdf = [filename_pdf '_chan' num2str(n) '_' num2str(snr) 'dB_snr.pdf']     % add info about snr to the pdf filename
    
    
    
    % plot and format the figure
    
    figure
    contour(t, f,  s_dB, [0:-1*snr_step:-1*snr])
    set(gca,'YScale','log')
    xlabel('Time (sec)')
    ylabel('Freq (Hz)')
    
    text(2, fs/4, num2str(fs))
    
    axis([0 max(t) 10 upper_frequency_limit])
    grid on
    c = colorbar;
    ylabel(c, 'SNR rel. peak level (dB)');
    
    ax = gca;
    
    ax.XTick = [0 10 20 30 40 50 60];
    ax.XTickLabel = {'0', '10', '20', '30', '40', '50', '60'};
    
    ax.YTick      = [10 100 1000 10000];
    ax.YTickLabel = {'10', '100', '1k', '10k'};
    ax.FontSize = 8;
    ax.YMinorGrid = 'on';
    ax.XMinorGrid = 'off';
    title(filename_pdf(1:end-4), 'Interpreter', 'none');
    
    set(gca, 'FontName','Helvetica')
    set(get(gca,'Title'),'FontName','Helvetica')
    set(get(gca,'XLabel'),'FontName','Helvetica')
    set(get(gca,'YLabel'),'FontName','Helvetica')
    set(gcf, 'Color', 'w');
    
    if savefile_on == 1
        set(gcf, 'PaperPosition', [0 0 21 14]);    %Position plot at left hand corner.
        set(gcf, 'PaperSize', [20 14]);              %Set the paper width and height.
        saveas(gcf, fullfile(pathname, filename_pdf), 'pdf') %Save figure
    end
    
end

end




%%

function output = click_generator(audio_channels, duration_in_seconds, level_dBFS, fs, bit_depth)

% generate a signal with a 1-sample click per second at a fixed level.
%
% output = click_generator(audio_channels, duration_in_seconds, level_dBFS, fs, bit_depth)
%
% audio_channels        Number of audio channels
% duration_in_seconds   total length of sweep in seconds
% level_dBFS            Level of click in dB FS
% fs                    Sampling rate in Hz
% bit_depth             Word length for TPDF dithering and quantisation
%
% geoff martin, Bang & Olufsen 2017

duration_in_samples = duration_in_seconds * fs;
click_level_lin = 10^(level_dBFS/20);

output = zeros(duration_in_samples, audio_channels);
output(1:fs:end, :) = click_level_lin;



end





%%

function output = cosine_sweep(audio_channels, f_start, f_stop, duration_in_seconds, level_dBFS, fs, bit_depth, zero_padding, dirac_on)

% generate a TPDF-dithered, quantised cosine wave with a swept frequency at a fixed level.
% The sweep can be preceded by a single-sample dirac for the purposes of time-alignment
%
% output = cosine_sweep(audio_channels, f_start, f_stop, duration_in_seconds, level_dBFS, fs, bit_depth, zero_padding, dirac_on)
%
% audio_channels        Number of audio channels
% f_start               Start frequency in Hz
% f_stop                End frequency in Hz
% duration_in_seconds   total length of sweep in seconds
% level_dBFS            Level of signal in dB FS
% fs                    Sampling rate in Hz
% bit_depth             Word length for TPDF dithering and quantisation
% zero_padding          The length of silence before and after the sweep, and
%                       before the dirac in seconds
% dirac_on              1 = include dirac, 0 = do not include dirac
%
% geoff martin, Bang & Olufsen 2017






% make some constants
level_lin = 10^(level_dBFS/20);

t = 0:1/fs:duration_in_seconds;

duration_in_samples = length(t);




% create the cosine sweep
signal = level_lin * chirp(t, f_start, duration_in_seconds, f_stop, 'logarithmic');
signal=signal';

% make the TPDF dither signal
dither = rand(duration_in_samples, 1) - rand(duration_in_samples, 1);

% add TPDF dither to the bit_depth of the signal
lin_gain = 2^(bit_depth-1) - 2;   % leaves 2 LSB's of headroom: 1 for two's complement and 1 LSB of dither
output_dithered = round(signal * lin_gain + dither); % quantise with dither.

output_dithered = output_dithered / 2^(bit_depth-1); %scale back down to �1 max

if dirac_on == 0
    output_dithered = [zeros(zero_padding*fs, 1); output_dithered; zeros(zero_padding*fs, 1)];
else
    output_dithered = [zeros(zero_padding*fs, 1); 1; zeros(zero_padding*fs, 1); output_dithered; zeros(zero_padding*fs, 1)];
end

% make the signal 2-channel, dual mono
output_dithered = repmat(output_dithered, 1, audio_channels);

output = output_dithered;


end




function latency_analysis(y, pathname, filename)

% rough preliminary analysis of latency difference between two channels
% of an audio file vs. time.
%
% latency_analysis(y, pathname, filename)
%
% y             audio input (wav file)
% pathname      path of audio file
% filename      name of audio file
%
% geoff martin, Bang & Olufsen 2017


savefile_on = 1;    % set to 0 to just plot the result without saving a pdf


filename_pdf = filename;    % match the filename of the .wav input and the .pdf output
filename_pdf(end-3:end) = [];   % strip the suffix from the audio file's name
filename_pdf = [filename_pdf '_latency_drift.pdf']     % add info to the pdf filename




end


function y = rawaudioread(filename, channels, bitdepth, encoding)

% output = rawaudioread(filename, audio_channels, bit_depth, encoding)
%
% filename
% audio_channels: number of audio channels. default = 2
% bitdepth: the word length of the audio samples in bits. default = 16
% encoding: encoding type
%   0: Signed Bigendian PCM
%   Other types to be added in the future


fileID = fopen(filename);
rawdata = fread(fileID, sprintf('bit%d', bitdepth));


% de-interleave
for n = 1:channels
    y(:,n) = rawdata(n:channels:end);
end


end


%% Get Filter for cos_sweep analysis


function h = getFilter(Number)
    global h_raisedco Hlp1 Hlp6 b2

    
    List = char(1:11); %Enumerate numbers of filters
    
    
    if ~isempty(find(List == Number))
        disp('There is no filter corresponding to this number')
    end
    
    switch Number
        case '1'
        h = h_sinc;  %sinc filter with Kaiser window
        case '2' %Elliptic direct-form
        h = Hlp1;     
        case '3'
        h = Hlp6;    
        case '4'
        h = b2;    
            
    end         
        
end

