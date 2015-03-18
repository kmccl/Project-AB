
%This program trims the zero padding off the front and back of the stimulus
%and plots the stimulus with the original. Awesome functions written by CWB
%and stolen by KSM on 1/13/15

[FileName,PathName]=uigetfile('*.wav','Select the Subject .wav file for calibration')
[data, fs] = audioread([PathName FileName])
[data_trim] = threshclipaudio(data, eps, 'begin&end');
[x, y, lags] = align_timeseries(data, data_trim, 'xcorr', 'fsx', 10000, 'fsy', 10000, 'pflag', 2);
size(data_trim)
size(data)
plot(data)
[filename, pathname] = uiputfile('*.wav', 'Choose a file name'); 
outname = fullfile(pathname, filename); 
audiowrite(outname, data_trim, 10000,'bitspersample', 32); 
