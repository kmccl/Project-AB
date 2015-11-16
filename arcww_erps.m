% arcww_erps.m
% Script to calculate and plot the ERPs for each condition and each subject

%subjects = {'1' '3' '4' '5' '6' '7' '9' '10'}; 
%subjects = {'1' '3' '6' '9' '10'}; % 7 still has saccades in data.
subjects = {'7'};

group_dir = '/auto/iduna/kbacker/ARCWW/group/erps/';

sem_erps = cell(length(subjects),1);
spa_erps = cell(length(subjects),1);
neu_erps = cell(length(subjects),1);

for s = 1:length(subjects)
    [si] = arcww_subject_info(subjects{s});
    
    in_dir = si.out_path; % Preproc directory
    %root_fn = [subjects{s},'_postica_latebase_100_in_r_b_a_s_'];
    %root_fn = [subjects{s},'_post-new-ica_in_r_b_a_s_'];
    root_fn = [subjects{s},'_post-ica_noeog_in_r_b_a_s_'];
    out_dir = ['/auto/iduna/kbacker/ARCWW/subjects/',subjects{s},'/erps/'];
    if ~exist(out_dir,'dir')
        mkdir(out_dir);
    end
    
    % Loop through each condition:
    for c = 1:length(si.rc_labels)
        
        % Load pre-processed dataset corresponding to the condition:
        [EEG] = pop_loadset('filename',[root_fn,si.rc_labels{c},'.set'],'filepath',in_dir);
        
        % Loop through each channel:
        erps = [];
        for ch = 1:EEG.nbchan
            m = mean(EEG.data(ch,:,:),3); % Average across epochs, 3rd dimension            
            erps = [erps; m];
        end
        
        % Low-Pass filter the averaged waveform in each channel:
        [filt_erps,a,b]=filt_lp(erps, si.filter.filt_length, si.filter.lp_cutoff, si.resamp_rate, 0);
        %close;
        
        % Save filtered erps as a matlab file:
        out_fn = [out_dir,subjects{s},'_',si.rc_labels{c},'.mat'];
        save(out_fn,'filt_erps');
        
        % Add data to appropriate cell array:
        if strcmpi(si.rc_labels{c},'Semantic')
           sem_erps{s} = filt_erps;
        elseif strcmpi(si.rc_labels{c},'Spatial')
           spa_erps{s} = filt_erps;
        elseif strcmpi(si.rc_labels{c},'Neutral')
           neu_erps{s} = filt_erps; 
        end
        
    end % for c
    
    % Plot Individual Topography:
    figure
   
    % Convert Samples to msec and make them relative to the R-C onset:
    samps = [1:size(sem_erps{s},2)];
    secs = samps/si.resamp_rate;
    secs2 = secs + si.epoch.timelim(1);
    msecs = secs2*1000;
    rmsecs = round(msecs);
    
    % Let's plot 5 timepoints for each condition:
    % 180 ms, 600 ms, 980 ms, 1400 ms, and 1800 ms:
    t = [180 600 980 1400 1800];
    
    %idx = [1220 1325 1420 1525 1625];
    sidx = 1;
    for i = 1:length(t)
        idx = find(rmsecs==t(i));
        
        % Neutral:
        subplot(length(t),3,sidx);
        topoplot(neu_erps{s}(:,idx),EEG.chanlocs,'maplimits',[-5 5],'style','map');
        title(sprintf('Subject %s: Neutral Retro-Cue\n Time: %s msec',subjects{s},num2str(msecs(idx))))
        axis tight
        colorbar
        sidx = sidx + 1;
        
        % Semantic:
        subplot(length(t),3,sidx);
        topoplot(sem_erps{s}(:,idx),EEG.chanlocs,'maplimits',[-5 5],'style','map');
        title(sprintf('Subject %s: Semantic Retro-Cue\n Time: %s msec',subjects{s},num2str(msecs(idx))))
        axis tight
        colorbar
        sidx = sidx + 1;
        
        % Spatial:
        subplot(length(t),3,sidx);
        topoplot(spa_erps{s}(:,idx),EEG.chanlocs,'maplimits',[-5 5],'style','map');
        title(sprintf('Subject %s: Spatial Retro-Cue\n Time: %s msec',subjects{s},num2str(msecs(idx))))
        axis tight
        colorbar
        sidx = sidx + 1;
        
    end % for i
    
end % for s

% Loop through each condition and get the group average ERP for each
% channel:
group_sem_erps = zeros(size(sem_erps{1}));
group_spa_erps = zeros(size(spa_erps{1}));
group_neu_erps = zeros(size(neu_erps{1}));

for i = 1:3 % Loop through each condition
    for j = 1:size(sem_erps{1},1) % Loop through each channel
        temp_ch = [];
        
        for k = 1:length(subjects)
            if i == 1 % Semantic
                temp_ch = [temp_ch; sem_erps{k}(j,:)];
            elseif i == 2 % Spatial
                temp_ch = [temp_ch; spa_erps{k}(j,:)];
            elseif i == 3 % Neutral
                temp_ch = [temp_ch; neu_erps{k}(j,:)];
            end % if
        end % for k
        
        m = mean(temp_ch);
        
        if i == 1 % Sem
            group_sem_erps(j,:) = m; 
        elseif i == 2 % Spa
            group_spa_erps(j,:) = m;
        elseif i == 3 % Neu
            group_neu_erps(j,:) = m;
        end
        
    end % for j
end % for i

% Convert Samples to msec and make them relative to the R-C onset:
samps = [1:size(group_sem_erps,2)];
secs = samps/si.resamp_rate;
secs2 = secs + si.epoch.timelim(1);
msecs = secs2*1000;

% Save data:
save([out_dir,'Group_SemanticNewICA.mat'], 'group_sem_erps', 'msecs');
save([out_dir,'Group_SpatialNewICA.mat'], 'group_spa_erps', 'msecs');
save([out_dir,'Group_NeutralNewICA.mat'], 'group_neu_erps', 'msecs');

% Plot ERP Data:
% Make Movie of Topography:
idx = [1175:5:1660]; % samples from 0 to 1.94 seconds.
z = figure;
subplot(1,3,1);
topoplot(group_neu_erps(:,idx(1)),EEG.chanlocs,'maplimits',[-3 3],'style','map');
title(sprintf('Neutral Retro-Cue\n Time: %s msec',num2str(msecs(idx(1)))))
axis tight
colorbar

subplot(1,3,2);
topoplot(group_sem_erps(:,idx(1)),EEG.chanlocs,'maplimits',[-3 3],'style','map');
title(sprintf('Semantic Retro-Cue\n Time: %s msec',num2str(msecs(idx(1)))))
axis tight
colorbar

subplot(1,3,3);
topoplot(group_spa_erps(:,idx(1)),EEG.chanlocs,'maplimits',[-3 3],'style','map');
title(sprintf('Spatial Retro-Cue\n Time: %s msec',num2str(msecs(idx(1)))))
axis tight
colorbar

set(gca,'nextplot','replacechildren');
aviobj = avifile('All_Conds_LargeNewICA.avi','fps',1,'quality',100);
for i = 1:length(idx)
    subplot(1,3,1);
    topoplot(group_neu_erps(:,idx(i)),EEG.chanlocs,'maplimits',[-3 3],'style','map');
    title(sprintf('Neutral Retro-Cue\n Time: %s msec',num2str(msecs(idx(i)))))
    
    
    subplot(1,3,2);
    topoplot(group_sem_erps(:,idx(i)),EEG.chanlocs,'maplimits',[-3 3],'style','map');
    title(sprintf('Semantic Retro-Cue\n Time: %s msec',num2str(msecs(idx(i)))))

    
    subplot(1,3,3);
    topoplot(group_spa_erps(:,idx(i)),EEG.chanlocs,'maplimits',[-3 3],'style','map');
    title(sprintf('Spatial Retro-Cue\n Time: %s msec',num2str(msecs(idx(i)))))

    
    F(i) = getframe(z);%,[left bottom width height]);
    aviobj = addframe(aviobj,F(i));
end
aviobj = close(aviobj);
movie(z,F,1,1)
