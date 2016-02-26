function [si] = ab_subject_info(sub)
% function [si] = ab_subject_info(sub)
%
% Contains the settings for eeg_preproc for ARCWW experiment.
%
% K. Backer, 17 February 2013

% Convert subject id to string if it's not already:
if ~isstr(sub)
    sub = num2str(sub);
end

% Get Filenames of Data
si.data_path = ['/Users/kmccl/Documents/DATA/subjects/',sub,filesep,'raw/']; % Directory where raw data is stored.
fn_prefix = ['AB0',sub];
fn_suffix = '.cnt';

if strcmp(sub,'11')
    si.fns = {'' '' '' '' 'ARCWW-11-4.cnt'};
else  
    fn_suffix2 = {'_MMBF7' '_Filtered'};
    temp = {};
    for f = 1:numel(fn_suffix2)
        temp{f} = [fn_prefix,fn_suffix2{f},fn_suffix];
    end % for f
    si.fns = temp;
end

% Where to save preprocessing outputs:
% Forward slanting separater is for unix and maybe Macs?
si.out_path = ['/Users/kmccl/Documents/DATA/subjects/',sub,filesep,'preproc/Sept2015/']; % Output of Pre-processing steps.
%si.tf_out_path = ['/home/kbacker/data/ARCWW/subjects/',sub,filesep,'tf/']; % Output for TF analysis
si.erp_out_path = ['/Users/kmccl/Documents/DATA/subjects/',sub,filesep,'erps/']; % Output for ERPs

% Make Directories if they don't exist:
if ~exist(si.out_path,'dir')
    mkdir(si.out_path);
end

% if ~exist(si.tf_out_path,'dir')
%     mkdir(si.tf_out_path);
% end

if ~exist(si.erp_out_path,'dir')
    mkdir(si.erp_out_path);
end

% Downsample to 250 Hz... keep getting out of memory error.
si.resamp_rate = 250;

% Remove EOG channels?
si.remove_eog = 1; %1 for yes, 0 for no.

% Study Pre-processing settings
% Jobs in order of pre-processing steps:
% im = import data, resample
% m = merge datasets, remove EOG
% e = epoch 
% ar = auto-reject bad trials
% i = ica first time
% v = visually inspect components epoch by epoch, reject noisy trials
% i = ica second time, then reject actual components
% in = interpolate bad channels over specified time ranges
% r = re-reference--Done after ICA because Cz is not independent of the
% other channels because it is the average reference.
% b = baseline
% a = artifact rejection (using threshold)
% s = sort trials based on condition and/or responses
% This ends the general pre-processing.


% Merge datasets before filtering to avoid multiple edge effects.
%if ~exist([si.out_path,sub,'_im_m_e_ar_i.set'],'file')
if ~exist([si.out_path,sub,'_im_m_e_i.set'],'file')
%if ~exist([si.out_path,sub,'_im_e_r_b_a.set'],'file')
    si.jobs = {'im' 'm' 'e' 'i'}; % Stop after first ICA to look at ICA components and choose
    % which to get rid of in "plot_ica_comps.m"
    %si.jobs = {'im' 'e' 'r' 'b' 'a'};% 'ar'}; % Took out merging and ICA.
%elseif ~exist([si.out_path,sub,'_im_m_e_ar_i_v_i.set'],'file') 
%elseif ~exist([si.out_path,sub,'_im_e_ar_i_v_i.set'],'file') 
%    si.jobs = {'i'}; % Stop after second ICA to inspect ICA components.
    % After 2nd ICA, run plot_ica_comps.m and list the eye/artifactual components
    % below.
%elseif ~exist([si.out_path,sub,'_im_m_e_icacorr.set'],'file') % Remove ICs:
elseif ~exist([si.out_path,sub,'_im_m_e_icacorr.set'],'file') % Remove ICs:
    si.jobs = {'i' 'r' 'b' 'a' 's'};
%elseif ~exist([si.out_path,sub,'_im_m_e_icacorr_in_r.set'],'file') % Finish general pre-processing:
%elseif ~exist([si.out_path,sub,'_im_e_icacorr_in_r_b_a.set'],'file') % Finish general pre-processing:
%    si.jobs = {'in' 'r' 'b' 'a'};
% elseif ~exist([si.out_path,sub,'_im_e_in_r_b_a.set'],'file') % Finish general pre-processing:
%     %si.jobs = {'in' 'e' 'r' 'b' 'a'};
%     si.jobs = {'in' 'r' 'b' 'a'};
    % These last steps will be done separately for ERPs and TF:
    % 'b' 'a' 's'; % Jobs post-ICA for ERP and TF analysis
% else
% si.jobs = {};% 's'};
% % else
% % si.jobs = {};% 's'}; These were commented out, not sure if I did it or KB
% % did it
end

% Filter Settings (actually, don't filter before ICA):
si.filter = struct( ...
     'type','lp', ... % hp = high-pass, lp = low-pass
     'hp_cutoff',0.1, ... % Cut off for hp filter
     'lp_cutoff',30,... % Cut off for lp filter, for ERPs
     'filt_length',255,...%500,...%1000,...
     'plotflag',1);

% Epoch and Cut settings:
temp = [10 20];
temp2 = {};
for i = 1:length(temp)
    temp2{i} = num2str(temp(i));
end
si.epoch.typerange = temp2;

% For now, time-lock to the Retro-Cue Onset.
% When epoching, add extra time at the end, because at least in prior
% versions of EEGLAB, it cuts off the last sample (so instead of 2000 ms,
% it would end at 1996 ms, assuming an FS of 250 Hz (4 ms per sample).
si.epoch.timelim  = [-0.3 .600]; % in seconds, Calculated based on various ITI and SOAs

% ERP Baseline Settings:
%si.baseline.timerange = [-3000 -2700]; % in msec, Relative to RC Onset,
%gives a huge shift in activity. must baseline right before event.
si.baseline.timerange = [-300 0]; % For ERPs.
%si.baseline.tf_timerange = [-4700 -4400]; % Otherwise, too many trials get rejected during artifact thresholding.

% After pre-processing, trim the epochs for ERP analysis:
%si.epoch.erp_timelim = [si.baseline.timerange(1)/1000 1.9];

% ICA Components to Reject:
% New pre-processing procedure:
% im_m_e_ar_i_v_i, then apply weights to im_m_e and THEN REJECT
% COMPONENTS!!!!!
% commented out 11/10 at 1109am
if strcmpi(sub,'393')
    si.icacomps = [1 7 9]; 
elseif strcmpi(sub,'395')
    si.icacomps = [1 3]; 
elseif strcmpi(sub,'396')
    si.icacomps = [1 4];  
elseif strcmpi(sub,'397')
    si.icacomps = [2]; 
elseif strcmpi(sub,'398')
    si.icacomps = [1 3 13]; 
elseif strcmpi(sub,'399')
    si.icacomps = [7 11]; 
elseif strcmpi(sub,'401')
    si.icacomps = [2 5 7 9]; 
 elseif strcmpi(sub,'402')
    si.icacomps = [5 7 14 16];
elseif strcmpi(sub,'403')
     si.icacomps = [1 2];
elseif strcmpi(sub,'405')
     si.icacomps = [1 3];
elseif strcmpi(sub,'406')
     si.icacomps = [1 3];
elseif strcmpi(sub,'407')
     si.icacomps = [1 3];
elseif strcmpi(sub,'408')
     si.icacomps = [1 4];
elseif strcmpi(sub,'417')
     si.icacomps = [1 2];
elseif strcmpi(sub,'418')
     si.icacomps = [1 3 4];
elseif strcmpi(sub,'421')
     si.icacomps = [1 4];
elseif strcmpi(sub,'423')
     si.icacomps = [1 2];
elseif strcmpi(sub,'427')
     si.icacomps = [1];
elseif strcmpi(sub,'429')
     si.icacomps = [1 4];
elseif strcmpi(sub,'430')
     si.icacomps = [1 4 30];
elseif strcmpi(sub,'431')
     si.icacomps = [3 8];
elseif strcmpi(sub,'434')
     si.icacomps = [1 4];
elseif strcmpi(sub,'409')
     si.icacomps = [1 2]; 
elseif strcmpi(sub,'410')
     si.icacomps = [1 3]; 
elseif strcmpi(sub,'421')
     si.icacomps = [1 4];  
elseif strcmpi(sub,'422')
     si.icacomps = [1 3];
elseif strcmpi(sub,'423')
     si.icacomps = [1 3];
elseif strcmpi(sub,'425')
     si.icacomps = [1 3];
elseif strcmpi(sub,'441')
     si.icacomps = [1 4];  
elseif strcmpi(sub,'445')
     si.icacomps = [1 7];
elseif strcmpi(sub,'446')
     si.icacomps = [1 4 7];
elseif strcmpi(sub,'437')
    si.icacomps = [1 2];
elseif strcmpi(sub,'449')
     si.icacomps = [2 3];
elseif strcmpi(sub,'443')
     si.icacomps = [1 2 3];
elseif strcmpi(sub,'438')
     si.icacomps = [1 4];  
elseif strcmpi(sub,'448')
     si.icacomps = [1 5];
elseif strcmpi(sub,'450')
     si.icacomps = [1 7];
elseif strcmpi(sub,'433')
     si.icacomps = [2];
elseif strcmpi(sub,'436')
     si.icacomps = [1 2];
elseif strcmpi(sub,'453')
     si.icacomps = [1 3];
elseif strcmpi(sub,'455')
     si.icacomps = [1 6];
elseif strcmpi(sub,'456')
     si.icacomps = [1 4];
elseif strcmpi(sub,'458')
     si.icacomps = [1 3];
end

% Threshold Artifact Rejection Settings:
% Need to do artifact thresholding 2x:
% Once for the post-retro-cue activity and once for the long baseline
% needed for time-frequency analysis.
si.thold1 = struct(...
    'typerej',1,... % 1 = raw data, 0 = ica
    'chans',[1:61],...% EOG Removed %'chans',[1:65],... % electrodes to examine
    'lowthresh',-150,... % lower t-hold in microvolts
    'upthresh',150,...% upper t-hold in microvolts
    'starttime',-.3,...% time to start in epoch (seconds)
    'endtime',0.6,... % time to end in epoch (seconds)
    'superpose',0,...
    'reject',1);

% This for TF, baseline artifact correction
% si.thold2 = struct(...
%     'typerej',1,... % 1 = raw data, 0 = ica    
%     'chans',[1:61],... %EOG Removed %'chans',[1:65],... % electrodes to examine
%     'lowthresh',-150,... % lower t-hold in microvolts
%     'upthresh',150,...% upper t-hold in microvolts
%     'starttime',-4.2,...% time to start in epoch (seconds)
%     'endtime',-2.7,...
%     'superpose',0,...
%     'reject',1);

% Trial Sorting Settings:
si.codes = {'10' '20'};
si.labels  = {'SPL' 'SL'};

% Interpolation Information:

% % Default epoch ranges:
% si.epochrange = {[1:144] [145:288] [289:432] [433:576]}; % If no trials
% %rejected before ICA.
% 
% if strcmpi(sub, '1')
%     si.epochrange = {[216:234] [235:345] [346:541] [542:573]};
%     si.badchans = {[27] [9 27] [9 27 22] [22]};
%     
% elseif strcmpi(sub, '3')
%     si.epochrange = {[1:576]};
%     si.badchans = {[9 27]};
%     
% elseif strcmpi(sub, '4')
%     si.epochrange = {[1:72] [73:288] [289:576]};
%     si.badchans = {[38 40 52] [38 40] [38]};
%     
% elseif strcmpi(sub, '5')
%     si.epochrange = {[1:72] [73:144] [145:287]};
%     si.badchans = {[14 18 21 27 28] [14 27 28] [14 28]};
%         
% elseif strcmpi(sub, '6')
%     si.epochrange = {[1:72] [73:211] [212:287] [288:360] [361:576]};
%     si.badchans = {[5 14 15 23 47 38] [14 23] [14 23] [14 27] [14]};
%     
% elseif strcmpi(sub, '7')
%     si.epochrange = {[1:71] [72:399] [400:547] [548:574]};
%     si.badchans = {[13 17 26 39 52 56] [39] [39 25] [39]};
% 
% elseif strcmpi(sub, '9')
%     si.epochrange = {[1:72] [72:200] [201:220] [221:290]};
%     si.badchans = {[7 26 42 8 16] [8 16] [8 16 26] [8]};
%     
% elseif strcmpi(sub, '10') 
%     si.epochrange = {[1:144]};
%     si.badchans = {[11]};
%     
% elseif strcmpi(sub,'11')
%     si.epochrange = {[1:72] [73:144]};
%     si.badchans = {[17 39] [17 39 54]};
%     
% elseif strcmpi(sub,'13')
%     si.epochrange = {};
%     si.badchans = {};
%     
% elseif strcmpi(sub,'14')
%     si.epochrange = {[73:288] [380:513]};
%     si.badchans = {[59] [2]};
%     
% elseif strcmpi(sub,'15')
%     si.epochrange= {[1:251] [252:390] [391:438] [439:575]};
%     si.badchans = {[17] [17 28] [17] [17 29]};
%     
% elseif strcmpi(sub,'16')
%     si.epochrange = {[144:250] [261:283] [284:567]};
%     si.badchans = {[1] [59] [1]};
%     
% elseif strcmpi(sub,'17')
%     si.epochrange = {[1:21] [22:34] [35:72] [73:216] [217:504]};
%     si.badchans = {[25] [24] [24 34] [34 23] [23 4]};
% 
% elseif strcmpi(sub,'8') % Noisy, trying to salvage
%     si.epochrange = {[1:72] [73:140] [141:170] [171:210] [211:281] [290:381] [382:462]};
%     si.badchans = {[15 16 23 42 45 54 60] [15 16 41 56 60] [15 60] [15 16 60] [16] [37] [10 27 37]};
% elseif strcmpi(sub,'12')
%     si.epochrange = {[1:72] [73:263] [264:345] [346:574]};
%     si.badchans = {[4 30 45 47 48 53] [30] [30 53] [30 53 17]};
% end % if 
    

