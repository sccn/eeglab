% function EEG = pop_rejcont(EEG);

function [EEG selectedregions com ] = pop_rejcont(EEG, varargin);

com = '';
if nargin < 1
    help pop_rejcont;
    return;
end;

if nargin < 2
    firstelec = 'EXG1';                % first non EEG channel
    % take all scalp electrodes
    % -------------------------
    if ~isempty(EEG.chanlocs)
        indelec = strmatch( firstelec, { EEG.chanlocs.labels });
        
        if isempty(indelec), elecrange = 1:EEG.nbchan;
        else                 elecrange = 1:(indelec-1);
        end;
    else
        elecrange = 1:EEG.nbchan;
    end;
    elecrange = deblank(vararg2str(elecrange));
    elecrange = elecrange(2:end-1);
    
    promptstr = { 'Channel range' ...
                  'Frequency range (Hz)' ...
                  'Frequency threshold in dB' ...
                  'Epoch segment length (s)' ...
                  'Minimum number of contiguous epochs' ...
                  'Add trails before and after regions (s)' ...
                  };
    initstr = { elecrange '35 128' '10' '0.5' '4' '0.25' };
    result = inputdlg2(promptstr, 'Reject portions of continuous data - pop_rejcont', 1, initstr);
    if length( result ) == 0 return; end;
    
    options = { 'elecrange'     str2num(result{1}) ...
                'freqlimit'     str2double(result{2}) ...
                'threshold'     str2double(result{3}) ...
                'epochlength'   str2double(result{4}) ...
                'contiguous'    str2double(result{5}) ...
                'addlength'     str2double(result{6}) };
else
    options = varargin;
end;

opt = finputcheck(options, { 'threshold'   'real'   []    10;
                             'elecrange'   'real'   []    [1:EEG.nbchan];
                             'freqlimit'   'real'   []    [35 128];
                             'contiguous'  'real'   []    4;
                             'addlength'   'real'   []    0.25;
                             'epochlength' 'real'   []    0.5 }, 'pop_rejcont');
     
%EEG.event = [];
grouplen  = opt.contiguous/2*opt.epochlength*EEG.srate+1; % maximum number of points for grouping regions
color     = [ 0 0.9 0];             % color of rejection window

% compute power spectrum
% ----------------------
NEWEEG = EEG;
% average reference 
% NEWEEG.data(opt.elecrange,:) = NEWEEG.data(opt.elecrange,:)-repmat(mean(NEWEEG.data(opt.elecrange,:),1), [length(opt.elecrange) 1]);

[TMPNEWEEG]= eeg_regepochs(NEWEEG, 0.25, [0 opt.epochlength], NaN);
%[TMPNEWEEG indices] = pop_rejspec(TMPNEWEEG, 1, [1:64], -100, 15, 30, 45, 0, 0);
%rejepoch = find(indices);
tmp   = fft(TMPNEWEEG.data, [], 2);
freqs = linspace(0, TMPNEWEEG.srate/2, size(tmp,2)/2);
freqs = freqs(2:end); % remove DC (match the output of PSD)
tmp   = tmp(:,2:size(tmp,2)/2,:);
warning('off', 'MATLAB:log:logOfZero');
tmpspec  = 10*log10(abs(tmp).^2);  
warning('on', 'MATLAB:log:logOfZero');
tmpspec  = tmpspec - repmat( mean(tmpspec,3), [1 1 TMPNEWEEG.trials]);
specdata = tmpspec;

% apply threshold to average of all electrodes
% --------------------------------------------
%[I1 Irej NS Erej] = eegthresh( mean(specdata(opt.elecrange, :, :), 1), size(specdata,2), 1:length(opt.elecrange), -100, 15, [freqs(1) freqs(end)], 30, 45);
[I1 rejepoch NS Erej] = eegthresh( mean(specdata(opt.elecrange, :, :), 1), size(specdata,2), 1, -100, opt.threshold, [freqs(1) freqs(end)], opt.freqlimit(1), opt.freqlimit(2));
fprintf('%d regions selected for rejection\n', length(rejepoch));

% build the winrej array for eegplot
% ----------------------------------
winrej   = [];
if ~isempty(find(cellfun(@isempty, { TMPNEWEEG.event.epoch }) == 1))
    error('Some events are not associated with any epoch');
end;
allepoch = [ TMPNEWEEG.event.epoch ];
for index = 1:length(rejepoch)
    eventepoch = find( rejepoch(index) == allepoch );
    if strcmpi(TMPNEWEEG.event(eventepoch(1)).type, 'X')
        urevent = TMPNEWEEG.event(eventepoch(1)).urevent;
        lat     = TMPNEWEEG.urevent(urevent).latency;
        %if length(
        winrej  = [ winrej; lat lat+opt.epochlength*NEWEEG.srate-1 color ]; %Erej(:,index)'];
    else
        error('Wrong type for epoch');
    end;
end;
winrej(:,6:6+length(opt.elecrange)-1) = 0;

% remove isolated regions and merge others
% ----------------------------------------
merged = 0;
isolated = 0;
for index = size(winrej,1):-1:1
    if size(winrej,1) >= index && winrej(index,2) - winrej(index,1) > grouplen, winrej(index,:) = []; isolated = isolated + 1;
    elseif index == 1 && size(winrej,1) > 1 && winrej(index+1,1) - winrej(index,2) > grouplen, winrej(index,:) = []; isolated = isolated + 1;
    elseif index == size(winrej,1) && size(winrej,1) > 1 && winrej(index,1) - winrej(index-1,2) > grouplen, winrej(index,:) = []; isolated = isolated + 1;
    elseif index > 1 && size(winrej,1) > 1 && index < size(winrej,1) && winrej(index+1,1) - winrej(index,2) > grouplen && ...
            winrej(index,1) - winrej(index-1,2) > grouplen
        winrej(index,:) = [];
        isolated = isolated + 1;
    elseif index < size(winrej,1) && size(winrej,1) > 1 && winrej(index+1,1) - winrej(index,2) <= grouplen
        winrej(index,2) = winrej(index+1,2);
        winrej(index+1,:) = [];
        merged = merged + 1;
    end;
end;
fprintf('%d regions merged\n', merged);
fprintf('%d regions removed\n', isolated);

% add time before and after each region
% -------------------------------------
for index = 1:size(winrej,1)
    winrej(index,1) = max(1,         winrej(index,1)-opt.addlength);
    winrej(index,2) = min(EEG.pnts,  winrej(index,2)+opt.addlength);
end;

% plot result
% -----------
selectedregions = winrej(:,1:2);
command = 'EEG = pop_select(EEG, ''nopoint'', TMPREJ(:,1:2)); [ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET, ''study'', ~isempty(STUDY)+0); eeglab redraw';
if nargin < 2
    eegplot(NEWEEG.data(opt.elecrange,:), 'srate', NEWEEG.srate, 'winrej', winrej, 'command', command, 'events', EEG.event);
end;
EEG = NEWEEG;
com = sprintf('EEG = pop_rejcont(EEG, %s);', vararg2str(options));
