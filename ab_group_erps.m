% ab_group_erps.m
% Script to calculate and plot the ERPs for each condition and each group

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Change the Subject numbers:
% subjects = {'1' '3' '4' '5' '6' '7' '9' '10'}; 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% conds = {'SPL' 'SL'}; % Don't change this order... it will screw things up when 
% % doing the group averages...
% 
% sl_erps = cell(length(subjects),1); % Cell to hold the SL ERPs
% spl_erps = cell(length(subjects),1); % Cell to hold the Fixed SPL ERPs
% 
% for s = 1:length(subjects)
%     [si] = ab_subject_info(subjects{s});
%     
%     in_dir = si.out_path; % Preproc directory
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Change this to where you want to save the ERPS for each subject:
%     % It makes the directory below.... no need to make it manually.
%     out_dir = ['/auto/iduna/kbacker/ARCWW/subjects/',subjects{s},'/erps/'];
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     if ~exist(out_dir,'dir')
%         mkdir(out_dir);
%     end
%     
%     % Loop through each condition:
%     for c = 1:length(conds)
%         
%         % Change the ending of the root_fn to how your files are named             
%         root_fn = [subjects{s},'_',conds{c},'_post-ica_noeog_in_r_b_a_s_'];
%         
%         % Load pre-processed dataset corresponding to the condition:
%         [EEG] = pop_loadset('filename',[root_fn,'.set'],'filepath',in_dir);
%         
%         % Loop through each channel:
%         erps = [];
%         for ch = 1:EEG.nbchan
%             m = mean(EEG.data(ch,:,:),3); % Average across epochs, 3rd dimension            
%             erps = [erps; m];
%         end
%         
%         % Low-Pass filter the averaged waveform in each channel:
%         [filt_erps,a,b]=filt_lp(erps, si.filter.filt_length, si.filter.lp_cutoff, si.resamp_rate, 0);
%         %close;
%         
%         % Save filtered erps as a matlab file:
%         out_fn = [out_dir,subjects{s},'_',si.rc_labels{c},'.mat'];
%         save(out_fn,'filt_erps');
%         
%         % Add data to appropriate cell array:
%         if strcmpi(conds{c},'SL')
%            sl_erps{s} = filt_erps;
%         elseif strcmpi(conds{c},'SPL')
%            spl_erps{s} = filt_erps;    
%         end        
%     end % for c
%     
% %     % Plot Individual Topography:
% %     figure
% %    
% %     % Convert Samples to msec and make them relative to the R-C onset:
% %     samps = [1:size(sem_erps{s},2)];
% %     secs = samps/si.resamp_rate;
% %     secs2 = secs + si.epoch.timelim(1);
% %     msecs = secs2*1000;
% %     rmsecs = round(msecs);
% %     
% %     % Let's plot 5 timepoints for each condition:
% %     % 180 ms, 600 ms, 980 ms, 1400 ms, and 1800 ms:
% %     t = [180 600 980 1400 1800];
% %     
% %     %idx = [1220 1325 1420 1525 1625];
% %     sidx = 1;
% %     for i = 1:length(t)
% %         idx = find(rmsecs==t(i));
% %         
% %         % Neutral:
% %         subplot(length(t),3,sidx);
% %         topoplot(neu_erps{s}(:,idx),EEG.chanlocs,'maplimits',[-5 5],'style','map');
% %         title(sprintf('Subject %s: Neutral Retro-Cue/n Time: %s msec',subjects{s},num2str(msecs(idx))))
% %         axis tight
% %         colorbar
% %         sidx = sidx + 1;
% %         
% %         % Semantic:
% %         subplot(length(t),3,sidx);
% %         topoplot(sem_erps{s}(:,idx),EEG.chanlocs,'maplimits',[-5 5],'style','map');
% %         title(sprintf('Subject %s: Semantic Retro-Cue/n Time: %s msec',subjects{s},num2str(msecs(idx))))
% %         axis tight
% %         colorbar
% %         sidx = sidx + 1;
% %         
% %         % Spatial:
% %         subplot(length(t),3,sidx);
% %         topoplot(spa_erps{s}(:,idx),EEG.chanlocs,'maplimits',[-5 5],'style','map');
% %         title(sprintf('Subject %s: Spatial Retro-Cue/n Time: %s msec',subjects{s},num2str(msecs(idx))))
% %         axis tight
% %         colorbar
% %         sidx = sidx + 1;
% %         
% %     end % for i
%     
% end % for s
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part does grand averaging:
%
% Loop through each condition and get the group average ERP for each
% channel:
groups = {'HLA' 'HLU' 'NH'}; % don't change this order...
subjects = {[398 401 403 405 430 431 434] [396 397 399 402 406 408 409 410 421 427 446 449 448] [393 395 407 417 418 422 423 425 429 441 445 437 443 438]};
conds = {'SPL' 'SL'}; % Don't change this order... it will screw things up when 
% doing the group averages...

% Setup main directory to load in data from (don't include the sub-directories for
% each group or condition.. those sub-directories are taken care of below.):
in_dir = '/Users/kmccl/Documents/DATA/subjects/Groups/';

% Load in some ERPs to get the size of the ERP matrix:
load([in_dir,'HLA/SL/SL398_.mat']);
% Whatever this variable that is loaded is called... replace wherever it
% says variable_name below...
hla_sl_erps = zeros(size(filt_erps));
hla_spl_erps = zeros(size(filt_erps));

hlu_sl_erps = zeros(size(filt_erps));
hlu_spl_erps = zeros(size(filt_erps));

nh_sl_erps = zeros(size(filt_erps));
nh_spl_erps = zeros(size(filt_erps));


% Loop through Each Group
for g = 1:length(groups)
    % Loop through each condition
    for c = 1:length(conds)        
        % load in data from all subjects in this group:
        sub_data = {};
        for s = 1:length(subjects{g})
            % The command below will load in the mat file for each
            % subject... whatever this variable that is loaded is called...
            % replace wherever it says variable_name2 below
            load([in_dir,groups{g},filesep,conds{c},filesep,...
                conds{c},num2str(subjects{g}(s)),'_.mat']);
            sub_data{s} = filt_erps;            
        end % for s
        
       % Loop through each channel:
        for j = 1:size(sub_data{1},1) 
            temp_ch = [];
            
            % Get channel j's ERP for each subject and put into temp_ch
            for s = 1:length(subjects{g})
                temp_ch = [temp_ch; sub_data{s}(j,:)];
            end % for k
            % Take the mean... aka grand average:
            m = mean(temp_ch);
            
            
            % Add the grand average to the right variable, depending on the
            % condition and group:
            if g == 1 % HLA
                if c == 1 % SPL
                    hla_spl_erps(j,:) = m;
                elseif c == 2 % SL
                    hla_sl_erps(j,:) = m;
                end
            elseif g == 2 % HLU
                 if c == 1 % SPL
                    hlu_spl_erps(j,:) = m;
                elseif c == 2 % SL
                    hlu_sl_erps(j,:) = m;
                end
            elseif g == 3 % NH
                 if c == 1 % SPL
                    nh_spl_erps(j,:) = m;
                elseif c == 2 % SL
                    nh_sl_erps(j,:) = m;
                 end               
            end
            
        end % for j
    end % for c
end % for g


% Convert Samples to msec and make them relative to the R-C onset:
[si] = ab_subject_info(subjects{1}(1));
samps = [1:size(sub_data{1},2)];
secs = samps/si.resamp_rate;
secs2 = secs + si.epoch.timelim(1);
msecs = secs2*1000;

% Save data: 
% To save all in the main Groups directory:
out_dir = in_dir;
save([out_dir,'Group_HLA_SPL.mat'], 'hla_spl_erps', 'msecs');
save([out_dir,'Group_HLA_SL.mat'], 'hla_sl_erps', 'msecs');
save([out_dir,'Group_HLU_SPL.mat'], 'hlu_spl_erps', 'msecs');
save([out_dir,'Group_HLU_SL.mat'], 'hlu_sl_erps', 'msecs');
save([out_dir,'Group_NH_SPL.mat'], 'nh_spl_erps', 'msecs');
save([out_dir,'Group_NH_SL.mat'], 'nh_sl_erps', 'msecs');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Plot ERP Data:
% % Make Movie of Topography:
% idx = [1175:5:1660]; % samples from 0 to 1.94 seconds.
% z = figure;
% subplot(1,3,1);
% topoplot(group_neu_erps(:,idx(1)),EEG.chanlocs,'maplimits',[-3 3],'style','map');
% title(sprintf('Neutral Retro-Cue/n Time: %s msec',num2str(msecs(idx(1)))))
% axis tight
% colorbar
% 
% subplot(1,3,2);
% topoplot(group_sem_erps(:,idx(1)),EEG.chanlocs,'maplimits',[-3 3],'style','map');
% title(sprintf('Semantic Retro-Cue/n Time: %s msec',num2str(msecs(idx(1)))))
% axis tight
% colorbar
% 
% subplot(1,3,3);
% topoplot(group_spa_erps(:,idx(1)),EEG.chanlocs,'maplimits',[-3 3],'style','map');
% title(sprintf('Spatial Retro-Cue/n Time: %s msec',num2str(msecs(idx(1)))))
% axis tight
% colorbar
% 
% set(gca,'nextplot','replacechildren');
% aviobj = avifile('All_Conds_LargeNewICA.avi','fps',1,'quality',100);
% for i = 1:length(idx)
%     subplot(1,3,1);
%     topoplot(group_neu_erps(:,idx(i)),EEG.chanlocs,'maplimits',[-3 3],'style','map');
%     title(sprintf('Neutral Retro-Cue/n Time: %s msec',num2str(msecs(idx(i)))))
%     
%     
%     subplot(1,3,2);
%     topoplot(group_sem_erps(:,idx(i)),EEG.chanlocs,'maplimits',[-3 3],'style','map');
%     title(sprintf('Semantic Retro-Cue/n Time: %s msec',num2str(msecs(idx(i)))))
% 
%     
%     subplot(1,3,3);
%     topoplot(group_spa_erps(:,idx(i)),EEG.chanlocs,'maplimits',[-3 3],'style','map');
%     title(sprintf('Spatial Retro-Cue/n Time: %s msec',num2str(msecs(idx(i)))))
% 
%     
%     F(i) = getframe(z);%,[left bottom width height]);
%     aviobj = addframe(aviobj,F(i));
% end
% aviobj = close(aviobj);
% movie(z,F,1,1)
