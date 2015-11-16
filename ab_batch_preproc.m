% ab_batch_preproc.m
%
% Batch code to loop through the pre-processing for multiple subjects
% overnight. :)  work while you sleep.
%
% KSM & KCB. One late night in the lab.  Robby don't be jealous.

subjects = {'425' '441' '445' '446'};
%'441''445' '446''409' '410' '421' '422' 
%437-something wrong. Error is:
%Undefined function or variable "idx".
%Error in ab_eeg_preproc (line 168)
   %[EEG] = pop_select(EEG,'notrial',[idx:idx+9]); % take out the first 10 20's
%Error in ab_batch_preproc (line 14)
    %ab_eeg_preproc(subjects(s));
% %'409' '421' '422' '425'
%{'393' '395' '396' '397' '398' '399' '401' '402' '403' '405'...
%     '406' '407' '408' '409' '410' '417' '418' '421' '422' '423' '425' '427' ...
%     '429' '430' '431' '434'};
for s = 1:length(subjects)
    ab_eeg_preproc(subjects(s));   
end % for s

%