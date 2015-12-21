% This is an attempt to hack KBs ERP code for Project AB
% KM June 2015

subjects = {'453'};
group_dir = '/Users/kmccl/Documents/ProjectAB/group/erps';
ab_erps = cell(length(subjects),1);

for s = 1:length(subjects)
    [si] = ab_subject_info(subjects{s});
    in_dir = si.out_path; % Preproc directory
            root_fn = [subjects{s},'_im_m_e_icacorr_r_b_a_s_SL'];
             %root_fn = [subjects{s},'_im_m_e_icacorr_r_b_a_s_SPL'];
             out_dir = ['/Users/kmccl/Documents/DATA/subjects/',subjects{s},'/erps/SL'];
              %out_dir = ['/Users/kmccl/Documents/DATA/subjects/',subjects{s},'/erps/SPL'];
    if ~exist(out_dir,'dir')
        mkdir(out_dir);
    end
    % Load pre-processed dataset corresponding to the condition:
    % [EEG] = pop_loadset('filename',[root_fn,subjects{s},'.set'],'filepath',in_dir);
    [EEG] = pop_loadset('filename',[root_fn,'.set'],'filepath',in_dir);
    % [EEG] = pop_loadset('filename',root_fn,'395_SL_im_e_ar.set','filepath',in_dir);
    
    % Loop through each channel:
    erps = [];
    for ch = 1:EEG.nbchan
        m = mean(EEG.data(ch,:,:),3); % Average across epochs, 3rd dimension
        erps = [erps; m];
    end
    % Low-Pass filter the averaged waveform in each channel:
    [filt_erps,a,b]=filt_lp(double(erps), si.filter.filt_length, si.filter.lp_cutoff, si.resamp_rate, 0);
    % close;
    
    % Save filtered erps as a matlab file:
      out_fn = [out_dir,subjects{s},'_','.mat'];
    %out_fn = [out_dir,subjects{s},'_','.set'];
   
end


%save(out_fn,'erps');
save(out_fn,'filt_erps')

% Convert Samples to msec and make them relative to the R-C onset:
samps = [1:size(erps,2)];
secs = samps/si.resamp_rate;
secs2 = secs + si.epoch.timelim(1);
msecs = secs2*1000;
rmsecs = round(msecs);

% out_dir = in_dir;
% save([out_dir,'393_SPL.mat'],'filt_erps','msecs');

%Plot Individual responses
%figure,plot(rmsecs,erps(channel,:))