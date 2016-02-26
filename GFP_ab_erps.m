% load('/Users/kmccl/Documents/DATA/subjects/Groups/Group_HLA_SL.mat')
% load('/Users/kmccl/Documents/DATA/subjects/Groups/Group_HLA_SPL.mat')
% load('/Users/kmccl/Documents/DATA/subjects/Groups/Group_HLU_SL.mat')
% load('/Users/kmccl/Documents/DATA/subjects/Groups/Group_HLU_SPL.mat')
% load('/Users/kmccl/Documents/DATA/subjects/Groups/Group_NH_SL.mat')
% load('/Users/kmccl/Documents/DATA/subjects/Groups/Group_NH_SPL.mat')

% GFP_hla_sl_erps=var(hla_sl_erps);
% GFP_hla_spl_erps=var(hla_spl_erps);
% GFP_hlu_sl_erps=var(hlu_sl_erps);
% GFP_hlu_spl_erps=var(hlu_spl_erps);
% GFP_nh_sl_erps=var(nh_sl_erps);
% GFP_nh_spl_erps=var(nh_spl_erps);


dataDir = [pwd '/'];
group = {'NH','HLU','HLA'};
subgroup = {'SPL','SL'};
grandDir = [dataDir 'Groups/']

for g = 1:length(group)
    for s = 1:length(subgroup)
        subDir = [dataDir group{g} '/' subgroup{s} '/'];
        files = dir(subDir);
        fnames = {files.name};
        matFiles = fnames(~cellfun(@isempty,regexp(fnames,'_.mat')));
        for f = 1:length(matFiles)
            load([subDir  matFiles{f}])
            GFP_filt_erps = var(filt_erps)%calc GFP
            save([subDir  matFiles{f}(1:end-4) 'GFP.mat'],'GFP_filt_erps')
        end
    end
end
