%from the EEGlablist
%Stefan Debener
% if you have an epoched EEGLAB dataset, you get the GFP from using the ML
% std() of the averaged signal. This doesn't work for my project because
% the averaged data is in .mat files and not .set files :(
addpath(genpath('/Users/kmccl/Documents/MATLAB'));

conds = {'SL' 'SPL'}; 
groups = {'NH' 'HLU' 'HLA'}; % don't change this order...
subjects = {[393 395 407 417 418 423 429] [396 397 399 402 406 408 427] [398 401 403 405 430 431 434]}; % Add in the additional subjects you've tested.
comps = {'P1' 'N1' 'P2'}; 

% Setup main directory to load in Group data from:
in_dir = '/Users/kmccl/Documents/DATA/subjects/393/preproc/Sept2015/'
EEG = pop_loadset([in_dir,'393_im_m_e_icacorr_r_b_a_s_SL.set']); % need channel locations to plot topographies.

if length(size(EEG.data))==3
    gfp = std(mean(EEG.data,3));
end

figure;
plot(EEG.times,mean(EEG.data,3),'k');
hold on;
plot(EEG.times,std(mean(EEG.data,3)),'r', 'linewidth',3);


