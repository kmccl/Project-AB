% Script to load the erpdata saved from EEGLAB and do permutation tests
% directly on that... keep getting out of memory error when using the STUDY
% function calls.  super annoying.
%
% Revised on June 22, 2014, because somewhere the code became incompatible
% with the saved data format in the group ERP directory.

group_dir = '/auto/iduna/kbacker/ARCWW/group/erps/n=16/';
cntrsts = {{'Semantic-Inf' 'Spatial-Inf' 'Neutral'}};
%cntrsts = {{'Semantic-Inf' 'Spatial-Inf' 'Neutral'} {'Semantic-Neu' 'Spatial-Neu'}...
%    {'Left' 'Right'} {'Salient' 'Not Salient'}};
%dcntrsts = {{'Semantic-Inf - Semantic-Neu' 'Spatial-Inf - Spatial-Neu'}};

% Set path to EEGLAB:
addpath(genpath('/home/kbacker/data/home/EEGLAB/eeglab11_0_4_4b/'));
num_perms = 2000;
thold = 0.005; % uncorrected.

% Load in an EEG dataset with the channel locations... needed below for
% topoplot:
[si] = new_arcww_subject_info('1');
EEG = pop_loadset('filename', '1_im_m_e_icacorr_in_r150_b.set','filepath',si.out_path);

for c = 1:length(cntrsts)
    load([group_dir,'Contrast',num2str(c),'.mat']);
    
    % KB. Revisions:
    % This loads sub_data, which is a cell array (1, num_conds)
    % Within each condition cell, there is a separate cell for each
    % subject.  Need to get it into the format that statcond will accept:
    erpdata = cell(size(sub_data));
    num_ch = size(sub_data{1}{1},1);
    num_t = size(sub_data{1}{1},2);
    num_s = length(sub_data{1});
    for x = 1:length(erpdata) % Loop through each cell
        erpdata{x} = zeros(num_ch, num_t, num_s);
        for y = 1:length(sub_data{x})
            erpdata{x}(:,:,y) = sub_data{x}{y};            
        end % for y
    end % for x
    clear sub_data;
    
    [stats,df,pvals] = statcond(erpdata,'paired','on','method','perm','naccu',num_perms);
    
    % Now, plot ERP data with the statistically thresholded results
    % (p-values): 
    
    electrodechart61;
    figure
    
    % Loop through each channel:
    %for ch = 1:size(erpdata{1},2)
    for ch = 1:size(erpdata{1},1) % channel is 1   
        % compute means: erpdata = timepoints (550) x channels (61) x subjects
        % (14)       
        % Now, erpdata = channels (61) x timepoints (550) x subjects (16)
        for e = 1:length(erpdata)
            %grp_data{e} = mean(erpdata{e}(:,plot_idx{ch}{3},:),3);
            grp_data{e} = mean(erpdata{e}(plot_idx{ch}{3},:,:),3);
        end % for e
        
        % pvals is timepoints (550) x channels (61):
        % Now, pvals is channels x timepoints.
        % Create a mask:
        %chp = pvals(:,plot_idx{ch}{3});
        chp = pvals(plot_idx{ch}{3},:);
        
        % Create erptimes variable:
        
        pmask = erptimes(find(chp<thold));
        z = zeros(length(erptimes),1);
        f = find(chp<thold);
        z(f) = 1;
        subplot(num_rows, num_cols, plot_idx{ch}{2});
        
        % Plot Significant Timepoints:
        for k = 1:length(z)
            if z(k)==1 % Significant
                plot([erptimes(k) erptimes(k)],[-3.49 3.49],'y','LineWidth',3);
                hold on;
            end
        end 
        
        if length(grp_data)==3
            i_hdl = plot(erptimes,grp_data{2},'r','LineWidth',3); % Plot first condition
            hold on;
            i_hdl = plot(erptimes,grp_data{3},'b','LineWidth',3); % Plot Second condition
            
            i_hdl = plot(erptimes,grp_data{1},'k','LineWidth',3); % Plot third condition if applicable.
        end
        
        % RC onset line:
        plot([0 0], [-3.5 3.5],':k');
        
       
        set(i_hdl,'Xdata',erptimes);
        set(i_hdl,'ButtonDownFcn','copyaxis');
        
        set(gca,'XLim',[erptimes(1) erptimes(end)],'YLim',[-3.5 3.5],'TickDir','out');
        axis([erptimes(1) erptimes(end) -3.5 3.5])
        title(plot_idx{ch}{1});
        
        %grid on
    end % for ch
    %leg_labels = [cntrsts{c} 'Cue Onset' ['p < ',num2str(thold)]];
    %legend(leg_labels);
    
    % Make topoplot at specified times:
    %times_to_plot = [172 228 600 700];
    times_to_plot = [540 1460];
    %clims = {[-6 6] [-3 3] [-2 2] [-2 2]};
    clims = {[-2 2] [-2 2]};
    all_grp_data = cell(size(erpdata));
    % Get Group mean ERPs WITH CHANNELS IN THE CORRECT ORDER!!!
    for e = 1:length(erpdata)
        for ch = 1:size(erpdata{e},2)
            all_grp_data{e} = [all_grp_data{e}; mean(erpdata{e}(:,ch,:),3)'];
        end % for ch
    end % for e
            
    for t = 1:length(times_to_plot)
        idx = find(erptimes==times_to_plot(t));
        figure
        
        subplot(1,length(all_grp_data),1);
        topoplot(all_grp_data{2}(:,idx),EEG.chanlocs,'maplimits',clims{t},'style','map');
        title(sprintf('%s Retro-Cue\n Time: %s msec',cntrsts{c}{1},num2str(times_to_plot(t))))
        colorbar
        
        subplot(1,length(all_grp_data),2);
        topoplot(all_grp_data{3}(:,idx),EEG.chanlocs,'maplimits',clims{t},'style','map');
        title(sprintf('%s Retro-Cue\n Time: %s msec',cntrsts{c}{2},num2str(times_to_plot(t))))
        colorbar
        
        if length(all_grp_data)==3
            subplot(1,3,3);
            topoplot(all_grp_data{1}(:,idx),EEG.chanlocs,'maplimits',clims{t},'style','map');
            title(sprintf('%s Retro-Cue\n Time: %s msec',cntrsts{c}{3},num2str(times_to_plot(t))))
            colorbar
        end
        
    end % for t
    
end % for c




for c = 1:length(dcntrsts)
    load([group_dir,'DContrast',num2str(c),'_ERPs.mat']);
    %[stats,df,pvals] = statcond(erpddata,'paired','on','method','perm','naccu',num_perms);
    
    % Now, plot ERP data with the statistically thresholded results
    % (p-values): 
    
    electrodechart61;
    figure
    all_grp_ddata = cell(size(erpddata));
    % Loop through each channel:
    for ch = 1:size(erpddata{1},2)
        
        % compute means: erpdata = timepoints (550) x channels (61) x subjects
        % (14)
        
        for e = 1:length(erpddata)
            grp_ddata{e} = mean(erpddata{e}(:,plot_idx{ch}{3},:),3);
            all_grp_ddata{e} = [all_grp_ddata{e}; grp_ddata{e}'];
        end % for e
        
        % pvals is timepoints (550) x channels (61):
        % Create a mask:
        chp = pvals(:,plot_idx{ch}{3});
        pmask = erptimes(find(chp<thold));
        z = zeros(length(erptimes),1);
        f = find(chp<thold);
        z(f) = 1;
        subplot(num_rows, num_cols, plot_idx{ch}{2});
        
        % Plot Significant Timepoints:
        for k = 1:length(z)
            if z(k)==1 % Significant
                plot([erptimes(k) erptimes(k)],[-3.49 3.49],'y','LineWidth',3);
                hold on;
            end
        end 
        
        if length(grp_ddata)==3
            i_hdl = plot(erptimes,grp_ddata{2},'r','LineWidth',3); % Plot first condition
            hold on;
            i_hdl = plot(erptimes,grp_ddata{3},'b','LineWidth',3); % Plot Second condition
            
            i_hdl = plot(erptimes,grp_ddata{1},'k','LineWidth',3); % Plot third condition if applicable.
        elseif length(grp_data)==2
            i_hdl = plot(erptimes,grp_ddata{1},'r','LineWidth',3); % Plot first condition
            hold on;
            i_hdl = plot(erptimes,grp_ddata{2},'b','LineWidth',3); % Plot Second condition
        end
        
        % RC onset line:
        plot([0 0], [-3.5 3.5],':k');
        
       
        set(i_hdl,'Xdata',erptimes);
        set(i_hdl,'ButtonDownFcn','copyaxis');
        
        set(gca,'XLim',[erptimes(1) erptimes(end)],'YLim',[-3.5 3.5],'TickDir','out');
        axis([erptimes(1) erptimes(end) -3.5 3.5])
        title(plot_idx{ch}{1});
        
        %grid on
    end % for ch
    %leg_labels = [cntrsts{c} 'Cue Onset' ['p < ',num2str(thold)]];
    %legend(leg_labels);
    
    
     % Make topoplot at specified times:
    times_to_plot = [540 644 1068 1460 1516 1864];
    clims = {[-2 2] [-2 2] [-2 2] [-2 2] [-2 2] [-2 2]};
    all_grp_ddata = cell(size(erpddata));
    % Get Group mean ERPs WITH CHANNELS IN THE CORRECT ORDER!!!
    for e = 1:length(erpddata)
        for ch = 1:size(erpddata{e},2)
            all_grp_ddata{e} = [all_grp_ddata{e}; mean(erpddata{e}(:,ch,:),3)'];
        end % for ch
    end % for e
            
    for t = 1:length(times_to_plot)
        idx = find(erptimes==times_to_plot(t));
        figure
        
        subplot(1,length(all_grp_ddata),1);
        topoplot(all_grp_ddata{1}(:,idx),EEG.chanlocs,'maplimits',clims{t},'style','map');
        title(sprintf('%s \n Time: %s msec',dcntrsts{c}{1},num2str(times_to_plot(t))))
        colorbar
        
        subplot(1,length(all_grp_ddata),2);
        topoplot(all_grp_ddata{2}(:,idx),EEG.chanlocs,'maplimits',clims{t},'style','map');
        title(sprintf('%s \n Time: %s msec',dcntrsts{c}{2},num2str(times_to_plot(t))))
        colorbar
        
      
        
    end % for t
    
    % Make Movie:
    % Make Movie of Topography:
    times_to_plot = [0:4:1900]; % samples from 0 to 1.94 seconds.
    idx = find(erptimes==times_to_plot(1));
    clim = [-2 2];
    z = figure;
    subplot(1,length(all_grp_ddata),1);
    topoplot(all_grp_ddata{1}(:,idx),EEG.chanlocs,'maplimits',clim,'style','map');
    title(sprintf('%s \n Time: %s msec',dcntrsts{c}{1},num2str(times_to_plot(1))))
    axis tight
    colorbar
    
    subplot(1,length(all_grp_ddata),2);
    topoplot(all_grp_ddata{2}(:,idx),EEG.chanlocs,'maplimits',clim,'style','map');
    title(sprintf('%s \n Time: %s msec',dcntrsts{c}{2},num2str(times_to_plot(1))))
    axis tight
    colorbar
    
    
    set(gca,'nextplot','replacechildren');
    aviobj = avifile(['DContrast',num2str(c),'_April2013.avi'],'fps',1,'quality',100);
    for i = 1:length(times_to_plot)
        idx = find(erptimes==times_to_plot(i));
        subplot(1,length(all_grp_ddata),1);
        topoplot(all_grp_ddata{1}(:,idx),EEG.chanlocs,'maplimits',clim,'style','map');
        title(sprintf('%s \n Time: %s msec',dcntrsts{c}{1},num2str(times_to_plot(i))))
        
        
        subplot(1,length(all_grp_ddata),2);
        topoplot(all_grp_ddata{2}(:,idx),EEG.chanlocs,'maplimits',clim,'style','map');
        title(sprintf('%s \n Time: %s msec',dcntrsts{c}{2},num2str(times_to_plot(i))))
        
       
        F(i) = getframe(z);%,[left bottom width height]);
        aviobj = addframe(aviobj,F(i));
    end
    aviobj = close(aviobj);
    movie(z,F,1,1)
    
    
end % for c