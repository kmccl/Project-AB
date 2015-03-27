% Comparing NAL-R to totally flattened filters:
% From 250 to 6000 Hz.
freqs = [250 500 1000 2000 3000 4000 6000];
thresholds = {[10 10 40 50 60 70 85] [25 30 30 45 50 65 65] [15 15 15 40 55 55 35] ...
    [5 10 5 10 15 10 20] [15 10 5 0 10 10 15] [15 20 20 35 40 50 45] [20 25 40 50 50 40 50]...
    [5 15 15 50 65 70 75] [15 15 20 30 55 60 65] [25 25 40 65 70 80 80] [15 15 15 30 50 55 55]}; % various audiogram thresholds
constants = [-17 -8 1 -1 -2 -2 -2]; 
retspl = [15.5 8.5 3.5 6.5 5.5 1.5 -1.5];% yost & killion RETSPL for ER-3As

% Loop through each audiogram:
for t = 1:length(thresholds)
    
    % Original Method:
    spl_t = thresholds{t} + retspl; % SPL Thresholds.
    adj_spl_t = spl_t - min(spl_t); % Setting the minimum to 0.
    
    % NAL-R Method:
    % Equation from Textbook of Hearing Aid Amplification
    % page 374.  Robert Sandlin, Editor. 2000.
    X = 0.05*(thresholds{t}(2) + thresholds{t}(3) + thresholds{t}(4));
    for f = 1:length(freqs) % find the gain for each frequency:
        if thresholds{t}(f) > 0
            gain_dBHL(f) = X + 0.31*(thresholds{t}(f))+constants(f);
        %else
        %    gain_dBHL(f) = 0;
        end
        % In Samira's code, she applied 0 gain for thresholds 25 or less.
        % If we don't do this constraint, maybe we could apply this NAL-R
        % to all subjects, including NH.  Or maybe that's a bad idea???
    end
    % Convert gain from dBHL to dB SPL:
    spl_gain = gain_dBHL + retspl;
    spl_gain0 = spl_gain - min(spl_gain); % Making the smallest value 0.
    
    figure,subplot(1,2,1)
    plot(freqs,spl_t,'k')
    title(['Audiogram #',num2str(t)])
    ylabel('Thresholds (dB SPL)')
    xlabel('Frequency')
    
    subplot(1,2,2)
    plot(freqs,adj_spl_t,'k')
    hold on, plot(freqs,spl_gain,'b')
    hold on, plot(freqs,spl_gain0,'r')
    legend('Flattened','NAL-R','NAL-R Set to 0')
    title(['Audiogram #',num2str(t)])
    ylabel('Gain (dB SPL)')
    xlabel('Frequency')
    axis([0 freqs(end) -30 100])
end