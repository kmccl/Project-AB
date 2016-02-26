%This is a little script to plot out my group averages so that I don't
%waste so much time retyping this out!
%KSM 11/2015
%Check the ELP file for channel numbers

% figure1=figure,plot(msecs,nh_sl_erps(61,:)', '-b', 'linewidth',3)
% hold on
% plot(msecs,hlu_sl_erps(61,:)','-r','linewidth',3)
% hold on
% plot(msecs,hla_sl_erps(61,:)','-k','linewidth',3)
% plottools
% 
% figure2=figure,plot(msecs,nh_spl_erps(61,:)', '-b', 'linewidth',3)
% hold on
% plot(msecs,hlu_spl_erps(61,:)','-r','linewidth',3)
% hold on
% plot(msecs,hla_spl_erps(61,:)','-k','linewidth',3)
% plottools


%For Laplacian figures
figure3=figure,plot(msecs,GFP_nh_sl_erps, '-b', 'linewidth',3)
hold on
plot(msecs,GFP_hlu_sl_erps,'-r','linewidth',3)
hold on
plot(msecs,GFP_hla_sl_erps,'-k','linewidth',3)
plottools

% figure4=figure,plot(msecs,GFP_nh_spl_erps, '-b', 'linewidth',3)
% hold on
% plot(msecs,GFP_hlu_spl_erps,'-r','linewidth',3)
% hold on
% plot(msecs,GFP_hla_spl_erps,'-k','linewidth',3)
% plottools