function [out,a,b]=filt_lp(in,len,lp,fs,plotflag)
% hp = small number, lp = large number
if ~exist('plotflag','var')
    plotflag = 0;
end
Nyq = fs/2;
%f = [0 hp/Nyq hp/Nyq lp/Nyq lp/Nyq 1]; %bandpass
%a = [0 0 1 1 0 0]; % bandpass
%f = [0 hp/Nyq hp/Nyq 1]; % Highpass
%a = [0 0 1 1]; % Highpass
f = [0 lp/Nyq lp/Nyq 1]; % Lowpass
a = [1 1 0 0]; % Lowpass
b = firls(len,f,a);

for i=1:size(in,1)
    % KB Added 5 July 2015: Zero-pad the ERPs before filtering... let's see
    % if this affects the filtering... The length of the ERPs (225 time
    % points) has to be at least 3x the filter order (aka filter length),
    % which is set in ab_subject_info... changed from 500 to 255... 
    erp = in(i,:);
    zero_erp = [zeros(1,500) erp zeros(1,500)];
    %out(i,:)=filtfilt(b,1,in(i,:));
    filt_erp = filtfilt(b, 1, zero_erp);
    % Now, take off the extra zeros from the filtered ERP:
    out(i,:) = filt_erp(501:end-500);
end

if plotflag == 1
    freqz(b,a,0:fs/2,fs);
end