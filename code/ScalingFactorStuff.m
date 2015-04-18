[FileName,PathName]=uigetfile('*.wav','Select the Subject .wav file for scaling');
[data, fs] = audioread([PathName FileName]);
figure,plot(data)
data_abs = abs(data);
scaling_factor = .9999/max(data_abs)
scaling_comparison=(scaling_factor)-[197.3050]
db_comparison=[db(scaling_factor)]-[db(197.3050)]
data_scaled=data.*scaling_factor;
hold on, plot([data data_scaled])
legend('Original Sound', 'Scaled Sound')



[filename, pathname] = uiputfile('*.wav', 'Choose a file name'); 
outname = fullfile(pathname, filename); 
audiowrite(outname, data_scaled, 10000,'bitspersample', 32);

