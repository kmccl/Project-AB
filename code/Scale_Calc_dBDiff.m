%This file was created to 
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

