%% Audio Sample Rate Converter from 44.1kHz to 48 kHz
clear all; clc;
edit AudioSRC_44p1_48_kHz

%% Change here to redesign SRC with desired SR

%SRC.OutputSampleRate = ;
disp('-----------------------');
disp('SRC from ... to ... kHz');
disp('-----------------------');
disp('Filter Informations');
disp('-------------------');
info(SRC);


%% Compare to single-stage
% load SRC_singleStage;
% filterBuilder(SRC_singleStage);

%% 