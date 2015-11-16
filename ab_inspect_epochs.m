% arcww_inspect_epochs.m
%
% script to examine epoched/ICA-rejected data to see which channels need to
% be interpolated for each participant.
%
% K Backer, March 26, 2013
%Modified June 2015 for collaboration with McClannahan & Tremblay

% subs = {'396'};

% for s = 1:length(subs)
    [si] = ab_subject_info(subs{s});
    
    % Loop through each cnt file of EEG data
    for f = 1:numel(si.fns)
        
        % Filename prefix:
        if isempty(strfind(si.fns{f},'ilter'))
            root_fn = ['396_SPL_im_e'];%took out ica step in preproc code for visual inspection
%             root_fn = [sub,'_SPL__im_e_icacorr'];
%             root_fn = [subs,'_SPL__im_e'];%took out 
        else
            root_fn = ['396_SL_im_e'];%took out ica step in preproc code for visual inspection
%             root_fn = [sub,'_SL__im_e_icacorr'];
%             root_fn = [subs,'_SL__im_e'];%took out 
        end
        
        %root_fn = [subs{s},'_im_e_icacorr'];
        EEG = pop_loadset('filename',[root_fn,'.set'],'filepath',si.out_path);
        
        eegplot(EEG.data,'srate',EEG.srate,'events',EEG.event)
        keyboard
        
    end % f
% end % for s