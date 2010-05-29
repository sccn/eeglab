disp('function is development to automatically reject portions of continuous data - Arnaud Delorme, May 2010');

return

% function EEG = pop_rejcont(EEG);
% 
% if nargin < 1
%     help pop_rejcont;
%     return;
% end;
% 
% f = finputcheck(

EEG.event = [];
firstelec = 'EXG1';                % first non EEG channel
epochlen  = 0.5;                   % epoch length
threshold = 10;                     % upper threshold in dB
freqrange = [35 128];              % frequency range
grouplen = 2*epochlen*EEG.srate+1; % maximum number of points for grouping regions
color    = [ 0 0.9 0];             % color of rejection window

% take all scalp electrodes
% -------------------------
indelec = strmatch( firstelec, { EEG.chanlocs.labels });
if isempty(indelec), indelec = EEG.nbchan; end;
elecrange = 1:(indelec-1);

% compute power spectrum
% ----------------------
NEWEEG = EEG;
NEWEEG.data(elecrange,:) = NEWEEG.data(elecrange,:)-repmat(mean(NEWEEG.data(elecrange,:),1), [length(elecrange) 1]);
[TMPNEWEEG]= eeg_regepochs(NEWEEG, 0.25, [0 epochlen], NaN);
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
%[I1 Irej NS Erej] = eegthresh( mean(specdata(elecrange, :, :), 1), size(specdata,2), 1:length(elecrange), -100, 15, [freqs(1) freqs(end)], 30, 45);
[I1 rejepoch NS Erej] = eegthresh( mean(specdata(elecrange, :, :), 1), size(specdata,2), 1, -100, threshold, [freqs(1) freqs(end)], freqrange(1), freqrange(2));
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
        winrej  = [ winrej; lat lat+epochlen*NEWEEG.srate-1 color ]; %Erej(:,index)'];
    else
        error('Wrong type for epoch');
    end;
end;
winrej(:,6:6+length(elecrange)-1) = 0;

% remove isolated regions and merge others
% ----------------------------------------
merged = 0;
isolated = 0;
for index = size(winrej,1):-1:1
    
    if index == 1 && winrej(index+1,1) - winrej(index,2) > grouplen, winrej(index,:) = [];
    elseif index == size(winrej,1) && winrej(index,1) - winrej(index-1,2) > grouplen, winrej(index,:) = [];
    elseif index > 1 && index < size(winrej,1) && winrej(index+1,1) - winrej(index,2) > grouplen && ...
            winrej(index,1) - winrej(index-1,2) > grouplen
        winrej(index,:) = [];
        isolated = isolated + 1;
    elseif index < size(winrej,1) && winrej(index+1,1) - winrej(index,2) <= grouplen
        winrej(index,2) = winrej(index+1,2);
        winrej(index+1,:) = [];
        merged = merged + 1;
    end;
end;
fprintf('%d regions merged\n', merged);
fprintf('%d isolated regions removed\n', isolated);

% add time before and after each region
% -------------------------------------
for index = 1:size(winrej,1)
    winrej(index,1) = max(1,         winrej(index,1)-epochlen*EEG.srate/2);
    winrej(index,2) = min(EEG.pnts,  winrej(index,2)+epochlen*EEG.srate/2);
end;

% plot result
% -----------
command = 'EEG = pop_select(EEG, ''nopoint'', TMPREJ(:,1:2)); [ALLEEG EEG CURRENTSET LASTCOM] = pop_newset(ALLEEG, EEG, CURRENTSET, ''study'', ~isempty(STUDY)+0); eeglab redraw';
eegplot(NEWEEG.data(elecrange,:), 'srate', NEWEEG.srate, 'winrej', winrej, 'command', command);
