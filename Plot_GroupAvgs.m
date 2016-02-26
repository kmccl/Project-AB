%This is a little script to plot out my group averages so that I don't
%waste so much time retyping this out!
%KSM 11/2015
%Check the ELP file for channel numbers

figure1=figure,plot(msecs,nh_sl_erps(41,:)', '-b', 'linewidth',3)
hold on
plot(msecs,hlu_sl_erps(41,:)','-r','linewidth',3)
hold on
plot(msecs,hla_sl_erps(41,:)','-k','linewidth',3)
plottools

figure2=figure,plot(msecs,nh_spl_erps(41,:)', '-b', 'linewidth',3)
hold on
plot(msecs,hlu_spl_erps(41,:)','-r','linewidth',3)
hold on
plot(msecs,hla_spl_erps(41,:)','-k','linewidth',3)
plottools


%For Laplacian figures
% figure3=figure,plot(msecs,nh_sl_laps(14,:)', '-b', 'linewidth',3)
% hold on
% plot(msecs,hlu_sl_laps(14,:)','-r','linewidth',3)
% hold on
% plot(msecs,hla_sl_laps(14,:)','-k','linewidth',3)
% plottools
% 
% figure4=figure,plot(msecs,nh_spl_laps(14,:)', '-b', 'linewidth',3)
% hold on
% plot(msecs,hlu_spl_laps(14,:)','-r','linewidth',3)
% hold on
% plot(msecs,hla_spl_laps(60,:)','-k','linewidth',3)
% plottools