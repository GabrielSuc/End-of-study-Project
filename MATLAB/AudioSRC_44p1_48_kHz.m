%% Initialization
Fs_In = 44.1e3;
Fs_Out = 48e3;
BW=22.05e3;

SRC = dsp.SampleRateConverter('Bandwidth',BW,'InputSampleRate',Fs_In,'OutputSampleRate',Fs_Out);
[L,M] = getRateChangeFactors(SRC);
SamplesPerFrame = 10*M;

AudioFileReader44p1 = dsp.AudioFileReader('Muddy_Waters.aiff','SamplesPerFrame',SamplesPerFrame);
AudioFileWriter48 = dsp.AudioFileWriter('Muddy_Waters_48k.wav','FileFormat','WAV','DataType','int32',...
    'SampleRate',Fs_Out);
AudioPlayer = audioDeviceWriter('SampleRate',Fs_Out);

SpectrumAnalyzer1 = dsp.SpectrumAnalyzer('SpectrumType','Power','Window','Hann','SampleRate',Fs_Out,...
    'PlotAsTwoSidedSpectrum',false,'Title','44.1kHz Audio Input','FrequencyScale','Linear','Position',[55 696 411 275]);
SpectrumAnalyzer2 = dsp.SpectrumAnalyzer('SpectrumType','Spectrogram','SampleRate',Fs_Out,...
    'PlotAsTwoSidedSpectrum',false,'Title','48kHz Audio Input','Position',[55 696 411 275]);
disp('------------------------------')
disp('Filter Info. and Cost Analysis')
disp('------------------------------')
info(SRC)
cost(SRC)

% PLOT


[y,Fs] = audioread('Muddy_Waters.aiff');

input = yc;
output = SRC(input);

plot(n_samp/Fs,y(1:44100),'*'); 
hold on 
plot(n_samp/(Fs*L/M),output(1:44100),'o'); 
hold off;
xlabel('Time (sec)');
ylabel('Signal value');
legend('Original','Resampled')
title('Resampling with DSP.SRC & Their filter')


%% Stream Processing Loop
tic;
while toc < 15
    % Run for 15 sec
    AudioIn = step(AudioFileReader44p1);
    % Sample conversion from 44p1 to 48 kHz
    AudioOut = step(SRC, AudioIn);
    % Play and display converted audio output
    step(AudioPlayer,AudioOut);         % Play resulting audio
    step(SpectrumAnalyzer1,AudioOut);   % Show touput power spectrum
    step(SpectrumAnalyzer2,AudioOut);   % Show output spectrogram
    step(AudioFileWriter48,AudioOut);   % Write to file   
end

%% Terminate 
release(SRC);
release(AudioPlayer);
release(AudioFileReader44p1);
release(AudioFileWriter48);