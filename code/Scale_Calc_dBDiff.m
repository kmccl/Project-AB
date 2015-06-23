%This file was created to find the maximum scaling factor for each
%individualized file, get the dB and compare that to the "perfect audio"
%scaling factor dB to find the amount of change for the attenuator. 
%KSM April 2015
clear all
close all

[FileName,PathName]=uigetfile('*.wav','Select the Subject .wav file for scaling');
[data, fs] = audioread([PathName FileName]);
figure,plot(data)
data_abs = abs(data);
scaling_factor = .9999/max(data_abs)
db_comparison =[[db(197.3050)] - [db(scaling_factor)]]
data_scaled=data.*scaling_factor;
hold on, plot([data data_scaled])
legend('Original Sound', 'Scaled Sound')
%Write new sound to file
[filename, pathname] = uiputfile('*.wav', 'Choose a file name'); 
outname = fullfile(pathname, filename); 
audiowrite(outname, data_scaled, 10000,'bitspersample', 32);

