
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
signal = level_lin * chirp(t, f_start, duration_in_seconds, f_stop, 'logarithmic',-90);
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