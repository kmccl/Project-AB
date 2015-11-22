% Script to look at ICA components and choose which ones to reject.
% Because for some reason when pop_selectcomp is called from a function, it
% breaks when you try to plot individual components.
%
% K. Backer, August 2, 2012
%393 already done separately
%subjects = [395 396 397 398 399 401 402 403 405 406 407 408 409 410 417 ...
    %417 418 421 422 423 425 427 429 430 431 434];% 3 4 5 6 7 9 10 11 13 14 15 16 17]395 409 410 421 422 425 437 449 443 438 448 ;
subjects = [450];
study = 'ab';
pause on
for s = 1:length(subjects)
    sub = num2str(subjects(s));
    [si] = feval([study,'_subject_info'],sub);
    
    % Loop through each cnt file of EEG data
    for f = 1:numel(si.fns)
        
%         % Filename prefix:
%         if isempty(strfind(si.fns{f},'ilter'))
%             root_fn = [sub,'_SPL_'];
%         else
%             root_fn = [sub,'_SL_'];
%         end
        root_fn = sub;
        
        % Load ICA set
        %EEG = pop_loadset('filename',[sub,'_im_m_e_ar_i_v_i.set'],'filepath',si.out_path);
        EEG = pop_loadset('filename',[root_fn,'_im_m_e_i.set'],'filepath',si.out_path);
        
        % Plot the components for inspection:
        pop_selectcomps(EEG,[1:60])
        %pop_selectcomps(EEG,[1:32])
        
        % User, input components for rejection:
        pause
        comps = input('Please enter the components you would like to reject (as a matrix): ');
        save(fullfile(si.out_path,[sub,'_Comps_to_Reject']),'comps'); % Saves as a .mat file
        
        % Subtracts the rejected ICA components from the data.
        % Also, plots the original and ICA-processed data on
        % the same plot for inspection and confirmation before
        % rejecting:
        if ~isempty(comps)
            for i = 1:length(comps) % Subtract 1 component at a time in this step...
                % Not actually saving the ICA-corrected data here... just
                % making sure that the subtraction of each component doesn't
                % distort brain activity
                sprintf('Displaying activity without Component %s',num2str(comps(i)))
                temp = pop_subcomp(EEG,comps(i),1);
            end
            
            % Now look at the data will all components rejected:
            display('Displaying activity without All Selected Components')
            temp = pop_subcomp(EEG,comps, 1);
        end
        
        %pause
    end % for f
end % for s