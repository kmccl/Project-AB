function [recon_snd] = individ_filter_stim_noscaling(snd, fs, audiogram, frequencies)
% function [recon_snd] = individ_filter_stim_noscaling(snd, fs, audiogram, frequencies)
% 
% This code "spectrally enhances" sounds, based on an individual's hearing
% profile.
% 
% INPUTS:
% 1) snd = the sound vector that you want to scale
% 2) fs = sound sampling rate
% 3) audiogram = the person's audiogram values (dB SPL) (from spreadsheet)
% 4) frequencies = the frequencies tested in the audiogram.
% The sound and audiogram files are not loaded in this version of the code.
% Thus, both of these files must be loaded into MATLAB by the code that
% calls this function.
%
% OUTPUT:
% 1) recon_snd
% This is the filtered, unscaled sound.  To make a long story short, scaling
% introduces various problems, depending on how the scaling is done.  So we
% removed the scaling altoghter and fixed some other bugs.
% 
% This code will also plot the spectrograms of the original and adjusted
% stimuli for comparison.
%
% Originally created by KM and CWB, October 2014.
% Modified by KCB and CWB, March 2015.

% Compute one-sided amplitude spectrum of original sound
NFFT = 2*length(snd);
fft_snd = fft(snd, NFFT)/NFFT;  % two-sided spectrum
amp = 2*abs(fft_snd(1:NFFT/2+1));

% Frequencies of amplitude spectrum
f = fs/2*linspace(0,1, NFFT/2+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CWB: 10 March 2015.
% Need to multiply amplitude by NFFT again for filtering process, lest we
% make the sound very small.
amp = amp .* NFFT; % Undo normalization.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subtract the minimum audiogram value to create the filter shape:
% audiogram = audiogram - min(audiogram);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KB SLIGHTLY MODIFIED INTERPOLATION PROCEDURE, PER TALKING WITH CHRIS
% BISHOP ON 4 MARCH 2015.
% RATHER THAN FINDING F>=FREQUENCIES(I) & F<FREQUENCIES(I+1), FIND
% THE FREQUENCY IN F THAT IS CLOSEST TO THE AUDIOGRAM
% FREQUENCIES... WHICH COULD BE BELOW THE AUDIOGRAM FREQUENCY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fidxs = [];
for x = 1:length(frequencies)
    % Lines from Chris:
    [sf, I] = sort(abs(f - frequencies(x))); 
    fidxs(x) = I(1); % This returns the indices
end % for x

% Interpolate our audiogram
interp_audiogram = [];
for i=1:numel(frequencies) - 1
    
    % Find the frequencies in this range
    % KB changed this to look for the closest frequency!    
    if i == numel(fidxs)-1 % If fitting the second-to-last and last point,
       % Also Include the point that is closest to the highest frequency
       mask = f>=f(fidxs(i)) & f<=f(fidxs(i+1));
    else
       mask = f>=f(fidxs(i)) & f<f(fidxs(i+1));
    end
    
    % Use a linear model to fit the two points
    [p] = polyfit([frequencies(i) frequencies(i+1)], [audiogram(i) audiogram(i+1)], 1); % original
    interp_audiogram = [interp_audiogram; [polyval(p, f(mask))]'];
    
end % for i=1:numel(audiogram)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF KB MODIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot audiogram and interpolated audiogram
figure, hold on
plot(frequencies, audiogram, 'rs--', 'linewidth', 2)

% Find range in amplitude spectrum to scale
% If change the masking procedure above, must also change it here!
% Also, include the last point:
mask = f>=f(fidxs(1)) & f<=f(fidxs(end));
plot(f(mask), interp_audiogram, 'ko');

% Prepend and append the first/last filter value to the interpolated
% audiogram.
interp_audiogram = [interp_audiogram(1).*ones(find(mask, 1, 'first')-1,1); ...
    interp_audiogram; interp_audiogram(end).*ones(find(fliplr(mask), 1, 'first') -1,1)];

% Scale magnitude (filter)
amp_filt = amp .* db2amp(interp_audiogram);

% Replace DC with the original DC. This should force the mean to be the
% same as the original file
amp_filt(1) = amp(1);

% Make figure to compare amplitude spectrum of original sound and the
% newly-filtered version (amp_filt), magnitude units
figure, hold on
plot(f, amp_filt, 'r');%, 'linewidth',2);
plot(f, amp, 'k');% 'linewidth', 2);
legend('Filtered Sound', 'Original Sound');

% Plot the amplitude spectra of the filtered and original sounds, this time
% in dB units:
figure, hold on
plot(f, db(amp_filt), 'r', 'linewidth',2);
plot(f, db(amp), 'k', 'linewidth', 2);
legend('Filtered Sound', 'Original Sound');
xlabel('Frequency (Hz)')
ylabel('Amplitude (dB)');

% Make two-sided spectrum
amp_filt_full = [amp_filt; flipud(amp_filt(2:end-1))]./2;

% Get the phase angle
ang = angle(fft_snd);

% make the fft
recon_fft = amp_filt_full.*cos(ang) + amp_filt_full.*sin(ang)*1i;
recon_snd= real(ifft(recon_fft));
recon_snd = recon_snd(1:length(snd));

% Windowing to remove large steps at beginning and end
%   Hanning ramp (5 ms)
[recon_snd] = fade(recon_snd, fs, true, true, @hann, 0.005);

% PLOT THE RECONSTRUCTED AND ORIGINAL TIME SERIES (SOUNDS)
figure
hold on
plot([recon_snd snd])
legend('Filtered Sound', 'Original Sound')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KB ADDED - PLOT THE SPECTROGRAMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the spectrograms of both sounds for comparison:
win_length = 256;
noverlap = win_length/2; %128;
nfft = win_length*2; % default = 256

% Spectrogram of Original Sound:
[S,F,T,P]=spectrogram(snd,win_length,noverlap,nfft,fs);
% Convert Power to correct units and plot
Pamp = abs(S);
%figure,imagesc(T,F,db(Pamp))
figure,imagesc(T,F,Pamp)
axis xy; axis tight;
xlabel('Time')
ylabel('Frequency')
title('Original Sound')
%caxis([0 10])
colorbar

% ALSO, PLOT THE SPECTROGRAM OF THE FILTERED SOUND:
% Spectrogram of Filtered Sound:
[S3,F3,T3,P3]=spectrogram(recon_snd,win_length,noverlap,nfft,fs);
% Convert Power to correct units and plot
Pamp3 = abs(S3);
%figure,imagesc(T3,F3,db(Pamp3))
figure,imagesc(T3,F3,Pamp3)
axis xy; axis tight;
xlabel('Time')
ylabel('Frequency')
title('Individualized Sound NOT SCALED')
%caxis([0 max(max(Pamp))]) % Plots same scale as previous graph.
colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF NEW ADDITION, FOR PLOTTING SPECTROGRAMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%