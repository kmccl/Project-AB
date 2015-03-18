%This is a program that modifies a sound stimulus according to auditory
%thresholds. Created by CWB and KSM on 10/20/14

clear all
close all

% Load audiogram from XLS file
[FileName,PathName]=uigetfile('*.xls','Select the subject Audiogram')
[data,txt,raw] = xlsread([PathName FileName])
%[data, txt, raw] = xlsread('HLU_1.xls'); 

% Transpose data
data = data';

% Plot audiogram
plot(data(:,1), data(:,6), 'ro--')
set(gca, 'Ydir', 'reverse');
xlabel('Frequency (Hz)');
ylabel('dB SPL');

% Load MMBF7.wav (/ba/)
%[ba, fs] = wavread('MMBF7.wav'); %USE FOR KATE'S COMPUTER
[ba, fs] = audioread('MMBF7.wav'); %USE FOR TESTING
% Plot ba
plot(ba); 

% Compute one-sided amplitude spectrum of /ba/
NFFT = 2*length(ba); 
fft_ba = fft(ba, NFFT)/NFFT;  % two-sided spectrum
amp = 2*abs(fft_ba(1:NFFT/2+1));

% Frequencies of amplitude spectrum
f = fs/2*linspace(0,1, NFFT/2+1);

% Plot amplitude spectrum
figure, hold on
plot(f, amp); 

% Compute filter response function from audiogram
% figure, hold on
audiogram = data(1:6, 6); 
frequencies = [data(1:6,1)'];

% Now find the cloest frequencies to these audiogram frequencies and use
% those

% Catch in case there
% audiogram = audiogram - max(audiogram); 
audiogram = audiogram - min(audiogram);
% audiogram = audiogram - mean(audiogram); 
% audiogram = [audiogram(1); audiogram]; 

% Interpolate our audiogram 
% audiogram_interp = [audiogram(1) audiogram audiogram(end)]; 
interp_frequencies = [];
interp_audiogram = [];
for i=1:numel(frequencies) - 1 
    
    % Find the frequencies in this range
    mask = f>=frequencies(i) & f<frequencies(i+1); 
    
    % Use a linear model to fit the two points
    [p] = polyfit([frequencies(i) frequencies(i+1)], [audiogram(i) audiogram(i+1)], 1);
    interp_audiogram = [interp_audiogram; [polyval(p, f(mask))]'];
    % Find the number of points between the frequencies 
%     nsamps = numel(find(f >= frequencies(i) & f < frequencies(i+1)));
    
    % interpolate frequencies
%     interp_frequencies = [interp_frequencies; [linspace(frequencies(i), frequencies(i+1), nsamps)]'];
    
    % interpolate audiogram
    % interpolate frequencies
%     interp_audiogram = [interp_audiogram; [linspace(audiogram(i), audiogram(i+1), nsamps)]'];
    
end % for i=1:numel(audiogram)

% Remove max again from interpolated audiogram
%   CWB is not sure why, but somewhere in the interpolation process, we're
%   changing the max of the audiogram by 0.0019 dB. So force this to be 0
%   (again)
interp_audiogram = interp_audiogram - min(interp_audiogram); 

% Plot audiogram and interpolated audiogram
figure, hold on
plot(frequencies, audiogram, 'rs--', 'linewidth', 2);


% Find range in amplitude spectrum to scale 
mask = f >= frequencies(1) & f < frequencies(end); 
plot(f(mask), interp_audiogram, 'ko');

% Prepend and append the first/last filter value to the interpolated
% audiogram.
interp_audiogram = [interp_audiogram(1).*ones(find(mask, 1, 'first')-1,1); interp_audiogram; interp_audiogram(end).*ones(find(fliplr(mask), 1, 'first') -1,1)];

% Scale magnitude (filter)
% amp_filt = amp;
amp_filt = amp .* db2amp(interp_audiogram); 

% Replace DC with the original DC. This should force the mean to be the
% same as the original file
amp_filt(1) = amp(1); 

% Make figure to compare amplitude spectrum of original sound (ba) and the
% newly-filtered version (amp_filt)
figure, hold on
plot(f, amp, 'k', 'linewidth', 2); 
plot(f, amp_filt, 'r'); 

% amp_filt(mask) = amp_filt(mask) .* db2amp(interp_audiogram);

% Flip amplitude spectrum
amp_filt_full = [amp_filt; flipud(amp_filt(2:end-1))]./2; 

% Get the phase angle
ang = angle(fft_ba); 

% make the fft
recon_fft = amp_filt_full.*cos(ang) + amp_filt_full.*sin(ang);

% Create an amplitude spectrum with the 0 dB frequencies matched. Makes it
% easier to compare the relative filter shape than using an RMS scaled
% version of the sound, which will match the MEAN dB across all
% frequencies. 
mask = interp_audiogram == interp_audiogram(end); 
plot(f, db(amp), 'k', 'linewidth', 2);
plot(f, db(amp_filt), 'r');
legend('Original /ba/', 'Filtered /ba/'); 
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)'); 

recon_ba= real(ifft(recon_fft));
recon_ba = recon_ba(1:length(ba));

% Windowing to remove large steps are beginning and end
%   Hanning ramp (5 ms)
[recon_ba] = fade(recon_ba, fs, true, true, @hann, 0.005); 

% RMS normalize ba, but only in the "0 dB" region of the filter. 

recon_ba_sc = recon_ba.*(rms(ba)./rms(recon_ba)); %USE FOR TESTING
%recon_ba_sc = recon_ba.* (max(abs(ba))./max(abs(recon_ba))); %USE FOR KATE'S COMPUTER

figure
hold on
plot([ba recon_ba_sc])

%creates the new .wav file
[filename, pathname] = uiputfile('*.wav', 'Choose a file name'); 
outname = fullfile(pathname, filename); 
% wavwrite(recon_ba_sc, fs, 32,outname); %USE FOR KATE'S COMPUTER
audiowrite(outname, recon_ba_sc, fs, 'bitspersample', 32); 
%audiowrite(recon_ba, fs, 32, 'HLU_1.wav'); %USE FOR TESTING

