function EEG = kb_eeg_interp(ORIEEG, bad_elec, bad_trials)
% function  EEG = kb_eeg_interp(ORIEEG, bad_elec, bad_trials)
% Based on eeglab's eeg_interp function.
%
% Does spherical interpolation only on epoched data.
% Added a bad_trials input to specify the trials during which to interpolate
% the bad electrodes, so that you can interpolate just a portion of the
% recording, without having to do it all.
%
% bad_trials is a matrix of trial (epoch) numbers during which to interpolate.
% 
% Kristina Backer, July 30, 2012

EEG = ORIEEG;

badchans  = bad_elec;
goodchans = setdiff(1:EEG.nbchan, badchans);
oldelocs  = EEG.chanlocs;
%EEG       = pop_select(EEG, 'nochannel', badchans);
% KB -- Added trial input to select only the epochs over which to
% interpolate:
EEG = pop_select(EEG, 'nochannel', badchans, 'trial', bad_trials);

%  So that the data can be re-concatenated in the correct order:
if bad_trials(1) == 1 && bad_trials(end) == size(ORIEEG.data,3) % All interpolated
    EEG1 = ORIEEG;
elseif bad_trials(1) == 1 || bad_trials(end) == size(ORIEEG.data,3) % First or last section
    good_trials = setxor(bad_trials, [1:size(ORIEEG.data,3)]);
    EEG1 = pop_select(ORIEEG, 'trial', good_trials);
else % Middle Section:
    good_trials1 = [1:bad_trials(1)-1];
    good_trials2 = [bad_trials(end)+1:size(ORIEEG.data,3)];
    EEG1 = pop_select(ORIEEG, 'trial', good_trials1);
    EEG2 = pop_select(ORIEEG, 'trial', good_trials2);
end
    
    
EEG.chanlocs = oldelocs;

% find non-empty good channels
% ----------------------------
origoodchans = goodchans;
chanlocs     = EEG.chanlocs;
nonemptychans = find(~cellfun('isempty', { chanlocs.theta }));

[tmp indgood ] = intersect(goodchans, nonemptychans);
goodchans = goodchans( sort(indgood) );
datachans = getdatachans(goodchans,badchans);
badchans  = intersect(badchans, nonemptychans);
if isempty(badchans), return; end;

% Using spherical method for interpolation:
% get theta, rad of electrodes
% ----------------------------
tmpgoodlocs = EEG.chanlocs(goodchans);
xelec = [ tmpgoodlocs.X ];
yelec = [ tmpgoodlocs.Y ];
zelec = [ tmpgoodlocs.Z ];
rad = sqrt(xelec.^2+yelec.^2+zelec.^2);
xelec = xelec./rad;
yelec = yelec./rad;
zelec = zelec./rad;
tmpbadlocs = EEG.chanlocs(badchans);
xbad = [ tmpbadlocs.X ];
ybad = [ tmpbadlocs.Y ];
zbad = [ tmpbadlocs.Z ];
rad = sqrt(xbad.^2+ybad.^2+zbad.^2);
xbad = xbad./rad;
ybad = ybad./rad;
zbad = zbad./rad;

EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
[tmp1 tmp2 tmp3 badchansdata] = spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad, EEG.data(datachans,:));

EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);

tmpdata = zeros(length(bad_elec), EEG.pnts, EEG.trials);
tmpdata(origoodchans, :,:) = EEG.data;
tmpdata(badchans , :) = badchansdata;
EEG.data = tmpdata;
EEG.nbchan = size(EEG.data,1);
EEG = eeg_checkset(EEG);

% Merge EEG (interpolated section) and EEG1/EEG2 (good section(s)):
tmp = EEG;

% Re-merge to keep trial order:
if bad_trials(1) == 1 && bad_trials(end) == size(ORIEEG.data,3) % Interpolated all:
    [EEG] = tmp;
elseif bad_trials(1) == 1 % Interpolated first section:
    [EEG] = pop_mergeset(tmp, EEG1);
elseif bad_trials(end) == size(ORIEEG.data,3) % Interpolated last section:
    [EEG] = pop_mergeset(EEG1, tmp);
else % Interpolated middle section:
    [tmp1] = pop_mergeset(EEG1, tmp);
    [EEG] = pop_mergeset(tmp1, EEG2);
end
[EEG] = eeg_checkset(EEG);

% get data channels
% -----------------
function datachans = getdatachans(goodchans, badchans);
      datachans = goodchans;
      badchans  = sort(badchans);
      for index = length(badchans):-1:1
          datachans(find(datachans > badchans(index))) = datachans(find(datachans > badchans(index)))-1;
      end;
      
function [xbad, ybad, zbad, allres] = spheric_spline( xelec, yelec, zelec, xbad, ybad, zbad, values);

newchans = length(xbad);
numpoints = size(values,2);

%SPHERERES = 20;
%[x,y,z] = sphere(SPHERERES);
%x(1:(length(x)-1)/2,:) = []; xbad = [ x(:)'];
%y(1:(length(x)-1)/2,:) = []; ybad = [ y(:)'];
%z(1:(length(x)-1)/2,:) = []; zbad = [ z(:)'];

Gelec = computeg(xelec,yelec,zelec,xelec,yelec,zelec);
Gsph  = computeg(xbad,ybad,zbad,xelec,yelec,zelec);

% compute solution for parameters C
% ---------------------------------
meanvalues = mean(values); 
values = values - repmat(meanvalues, [size(values,1) 1]); % make mean zero

values = [values;zeros(1,numpoints)];
C = pinv([Gelec;ones(1,length(Gelec))]) * values;
clear values;
allres = zeros(newchans, numpoints);

% apply results
% -------------
for j = 1:size(Gsph,1)
    allres(j,:) = sum(C .* repmat(Gsph(j,:)', [1 size(C,2)]));        
end
allres = allres + repmat(meanvalues, [size(allres,1) 1]);


% compute G function
% ------------------
function g = computeg(x,y,z,xelec,yelec,zelec)

unitmat = ones(length(x(:)),length(xelec));
EI = unitmat - sqrt((repmat(x(:),1,length(xelec)) - repmat(xelec,length(x(:)),1)).^2 +... 
                (repmat(y(:),1,length(xelec)) - repmat(yelec,length(x(:)),1)).^2 +...
                (repmat(z(:),1,length(xelec)) - repmat(zelec,length(x(:)),1)).^2);

g = zeros(length(x(:)),length(xelec));
%dsafds
m = 4; % 3 is linear, 4 is best according to Perrin's curve
for n = 1:7
    if ismatlab
        L = legendre(n,EI);
    else % Octave legendre function cannot process 2-D matrices
        for icol = 1:size(EI,2)
            tmpL = legendre(n,EI(:,icol));
            if icol == 1, L = zeros([ size(tmpL) size(EI,2)]); end;
            L(:,:,icol) = tmpL;
        end;
    end;
    g = g + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
end
g = g/(4*pi); 