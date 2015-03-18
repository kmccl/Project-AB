% Modified for use on Project AB 
% KM 3/18/15

clear all 
close all

% Load audiogram from XLS file
[FileName,PathName]=uigetfile('*.xls','Select the subject Audiogram')
[data,txt,raw] = xlsread([PathName FileName])

% Transpose data:
data = data';
audiogram = data(1:7, 6); % 1:8 selects audiogram frequencies up to 8000Hz
% 5000 Hz--Adding in the mean of 4000Hz and 6000Hz for /ba/ that has a Nyquist of 5000 Hz
sprintf('Original min of audiogram: %s',num2str(min(audiogram)))
frequencies = [data(1:7,1)'];

% Plot audiogram
figure, plot(data(1:7,1), data(1:7,6), 'ro--');
set(gca, 'Ydir', 'reverse');
xlabel('Frequency (Hz)');
ylabel('dB SPL');

%call function that scales /ba/
[recon_snd] = individ_filter_stim_noscaling(snd, fs, audiogram, frequencies);
recon_snd_sc = recon_snd * 0.1;


 % See if amplitude at best freq in filtered sound correlates with
 % scale factors.
 %[sfbas,idx] = sort(fbas);
  %sscale_factors = scale_factors(idx);
%figure,scatter(sfbas,sscale_factors)
            
 % SAVE THE SOUND FILE OUTPUTTED:
 %[filename, pathname] = uiputfile('*.wav', 'Choose a file name');
%outname = fullfile(pathname, filename);
 %wavwrite(recon_snd_sc, fs, nbits,outname); %USE FOR KATE'S COMPUTER
 %audiowrite(outname, recon_snd_sc, fs, 'bitspersample', nbits);
%audiowrite(recon_ba, fs, 32, 'HLU_1.wav'); %USE FOR TESTING
            
 % Save Adjusted Sound in Subject Directory:
 outfn = [out_stim_dir2,d(dd).name];
wavwrite(recon_snd_sc, fs, nbits, outfn);
            
 % Close the figures that pop up:
close all
           