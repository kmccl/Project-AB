%This is an unnecessary program that reads in the stimulus I want, scales
%it, and writes it to a file KM 3/18/15

[ba, fs] = audioread('MMBF7.wav');
ba_scaled = ba * 0.045;

figure
hold on
plot([ba ba_scaled])
legend('Original Sound', 'Scaled Sound')

[filename, pathname] = uiputfile('*.wav', 'Choose a file name'); 
outname = fullfile(pathname, filename);
audiowrite(outname, ba_scaled, fs, 'bitspersample', 32);
