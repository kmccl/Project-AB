%function [all_rms] = sea_make_aud_stim(grp_id, sub_id, main_study_dir)
% function sea_make_aud_stim(sub_id)
% grp_id = a number (1 to 3)
% sub_id = a number (1 to 30)
%
% Code to make individualized stim for the SEA Project
% Loops through both real and practice stim.
% K. Backer, February & March 2015.

% Modified for use on Project AB 
% KM 3/18/15

%sub_dir = [main_study_dir,'Group',num2str(grp_id),filesep,num2str(sub_id),filesep];
%out_dirs = {'Filtered_Ba'};

% Load in Audiogram Spreadsheet for this subject:
%fn = [sub_dir,'Group',num2str(grp_id),'_',num2str(sub_id),'.xls'];
%[data, txt, raw] = xlsread(fn);

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

all_rms = [];
for o = 1:length(out_dirs)
    % Select appropriate in_stim_dir:
    if strcmpi(out_dirs{o},'Adjusted_Sounds')
        in_stim_dir = [main_study_dir,'Stimuli\Sounds\'];
    elseif strcmpi(out_dirs{o},'Practice_Adjusted_Sounds')
        in_stim_dir = [main_study_dir,'Stimuli\Practice\'];
    end
    
    out_stim_dir = [sub_dir,out_dirs{o},filesep];
    if ~exist(out_stim_dir,'dir')
        mkdir(out_stim_dir);
    end
    
    % Loop through each sound in Normed_Animal and Normed_Music directories:
    cats = {'Normed_Animal' 'Normed_Music'};
    %scale_factors = [];
    %fbas =  [];
    for c = 1:length(cats)
        
        % Get the sounds in this directory:
        d = dir([in_stim_dir,cats{c},filesep,'*.wav']);
        
        out_stim_dir2 = [out_stim_dir,cats{c},filesep];
        if ~exist(out_stim_dir2,'dir')
            mkdir(out_stim_dir2);
        end
        
        % Load each Auditory Stimulus:
        for dd = 1:length(d)
            snd_fn = [in_stim_dir,cats{c},filesep,d(dd).name];
            [snd, fs, nbits] = wavread(snd_fn);
            
            [recon_snd] = individ_filter_stim_noscaling(snd, fs, audiogram, frequencies);
            recon_snd_sc = recon_snd * 0.1;
            %scale_factors = [scale_factors; scale_factor];
            %fbas = [fbas; fba];
            
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
            
            all_rms = [all_rms; rms(recon_snd_sc)];
        end % for dd
    end % for t
end % for o