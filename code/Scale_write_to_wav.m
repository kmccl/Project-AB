%This is an unnecessary program that reads in the stimulus I want, scales
%it, and writes it to a file KM 3/18/15
%modified to mess with the new filtered files and figuring out the scaling
%and attenuator settings. KM 4/15/15

clear all
close all

[FileName,PathName]=uigetfile('*.wav','Select the Subject .wav file for calibration');
[data, fs] = audioread([PathName FileName]);

data_scaled = data.*197.3050;

max(abs(data_scaled))

figure
hold on
plot([data data_scaled])
legend('Original Sound', 'Scaled Sound')

[filename, pathname] = uiputfile('*.wav', 'Choose a file name'); 
outname = fullfile(pathname, filename); 
audiowrite(outname, data_scaled, 10000,'bitspersample', 32);


% [PerfAudioFilt, fs] = audioread('PerfectAudio_filtered_0pad.wav');
% PerfAudioFilt_scaled = PerfAudioFilt.*;
% 
% max(abs(PerfAudioFilt_scaled))
% 
% figure
% hold on
% plot([PerfAudioFilt PerfAudioFilt_scaled])
% legend('Original Sound', 'Scaled Sound')
% 
% [filename, pathname] = uiputfile('*.wav', 'Choose a file name'); 
% outname = fullfile(pathname, filename);
% audiowrite(outname, PerfAudioFilt_scaled, fs, 'bitspersample', 32);
