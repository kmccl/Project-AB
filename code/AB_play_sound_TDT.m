[FileName,PathName]=uigetfile('*.wav','Select the appropriate /ba/ stimulus')
% file2play = [PathName FileName];
circuit = fullfile(pwd, '..', 'RP_Files', 'AB.rcx'); 
ntrials = 400;
fslevel=2;  %  2.441406250000000e+04
trigdur=0.02; % in seconds
trigdelay=1.24 + 63; % Delay trigger by 1.24 ms.  (in msecs)
code = 10; 

% isi = [0 0];
isi = [2 2.25]; % in seconds--this specifies the jitter window
% [2 2.25] 

%% LOAD STIMULI
%   Load using AA_loaddata. Only accept wave/double/single data types. 
p.datatypes=[1 2]; % restrict to wav and single/double array
if exist('tfs', 'var') && ~isempty(tfs), p.fs=tfs; end % work around for test click
[time_series, fs]=audioread(file2play); 

%% INITIALIZE TDT
% ESTABLISH CONNECTION WITH RP 2.1
RP=actxcontrol('RPCo.x',[0 0 0 0]);
RP.ConnectRP2('USB', 1);

% CLEAR Control Object File (COF)
RP.ClearCOF; % First clear the COF

% LOAD COF AND SET SAMPLING RATE
RP.LoadCOFsf(circuit, fslevel); 

RP.Halt;

% Zero out Buffer
%   Prevents silly mistakes like playing left-over sound from previous
%   setup.
RP.WriteTagV('Snd', 0, zeros(1, RP.GetTagVal('WavSize'))); 

% GET SAMPLING RATE
FS = RP.GetSFreq;

%% RESAMPLE STIMULI
%   Resample to TDT sampling rate
time_series=resample4TDT(time_series, FS, fs);

% Flip dimensions for TDT's sake. 
%   TDT wants a 1 x N vector. 
time_series = time_series'; 

%% CLIPPING CHECK
%   Double check that no clipping has occurred after summation. If so, then
%   throw an error.
%       This means the normalization routine above failed
if max(abs(time_series))>1, error('Stimulus clipped'); end 

%% SET TDT TAG VALUES
%   Several of these must be converted to samples first. 
if ~RP.SetTagVal('WavSize', length(time_series)), error('Parameter not set!'); end
if ~RP.WriteTagV('Snd', 0, time_series), error('Data not sent to TDT!'); end
if ~RP.SetTagVal('TrigDur', round(trigdur*FS)), error('Parameter not set!'); end       % TRIGger DURation in samples
if ~RP.SetTagVal('trigdelay', trigdelay), error('Parameter not set!'); end       
if ~RP.SetTagVal('TrigCode', code), error('Parameter not set!'); end               % Trigger CODE

%% START THE CIRCUIT
RP.Run;

%% MAIN CONTROL LOOP
%   Main stimulus control loop.

% make an 'edit' uicontrol to show the sweeps
tdisplay = uicontrol('style','edit','string', sprintf('Trial Number\nTime to Completion'),'units','normalized','position',[.2 .6 .6 .3], 'Max', 3, 'Min', 1);

% % Create a pause function for Kate
% pauseKey = KbName('P');
% KbQueueCreate([]);
% while KbCheck; end % Wait until all keys are released.
% KbQueueStart([]);

for n=1:ntrials
    
    % Present stimulus/drop trigger
    RP.SoftTrg(1);
    tic;
    jit=rand(1)*(diff(isi))+isi(1); 
    JIT(n)=jit; 
    % Round estimated remaining time to tenths place
    set(tdisplay,'string',sprintf(['Trial Number: %d (of %d).\n' ...
    'Estimated Time to Completion: %.1d s'],n, ntrials, round((length(time_series)./FS + mean(isi))*(ntrials-n)*10)/10));
    
    % Loop until stimulus is finished
    while ~RP.GetTagVal('TrialStatus'), end 
    
    % Pause for a while
    pause(jit); 
    
    % Check keyboard
%     [ pressed, firstPress]=KbQueueCheck([]);
%     if firstPress(pauseKey)
%         input('Press enter to continue'); 
%     end %
    
    toc
end % for n=1:ntrials

% Close keyboard queue
% KbQueueRelease([])