function [smoothdata,filtwts] = ERPLAB_butter1(data,srate,locutoff,hicutoff,samples,filterorder)

%     function [smoothdata,filtwts] = ERPLAB_butter1(data,srate,locutoff,hicutoff,samples,filterorder)
%
%     data        - input dataset
%     srate       - sample frequency
%     locutoff    - lower edge of the frequency pass band (Hz)  {0 -> lowpass}
%     hicutoff    - higher edge of the frequency pass band (Hz) {0 -> highpass}
%     filterorder - length of the filter in points {default 3*fix(srate/locutoff)}
%     
%     Outputs:
%     EEGOUT      - output dataset
%     
%     Author: ERPLAB Team, Center for Mind & Brain
%             Universidad de California, Davis.

if nargin < 1
	help ERPLAB_butter1
	return
end

if exist('filtfilt','file') ~= 2
    disp('ERPLAB_butter1 error: cannot find the Signal Processing Toolbox');
    return
end

if isempty(data)
    disp('ERPLAB_butter1 error: cannot filter an empty dataset')
    return
end

if nargin < 4
    disp('ERPLAB_butter1 error: please, enter all arguments!')
    return
end

if locutoff == 0 && hicutoff == 0,
   error('Need both lower and higher edge for you want to do');
   %return
end

% numchan  = EEG.nbchan;        % number of channels
fnyquist = 0.5*srate;           % half sample rate 
% samples  = EEG.pnts;
trans    = 0.15;                % fractional width of transition zones
minfac   = 3;                   % this many (lo)cutoff-freq cycles in filter
[numchan frames] = size(data);


if locutoff > fnyquist
    error('Low cutoff frequency cannot be > srate/2');
    %return
end

if hicutoff > fnyquist
    error('High cutoff frequency cannot be > srate/2');
    %return
end

if filterorder*3 > frames          % filtfilt restriction
    fprintf('ERPLAB_butter1: filter order too high');
    error('samples must be at least 3 times the filter order.');
end

if locutoff == 0 && hicutoff > 0     % Low Pass Filter
    [b,a] = butter(filterorder,(hicutoff/fnyquist));
    
    %suggestedfiltorder = minfac*fix(srate/hicutoff);
    %fprintf('Suggested filter order  is %d. ',suggestedfiltorder);
    
    if ((1+trans)*hicutoff/fnyquist) > 1
        error('high cutoff frequency too close to Nyquist frequency');
    end
    
elseif locutoff >0 && hicutoff == 0     % High Pass Filter
    [b,a] = butter(filterorder,(locutoff/fnyquist),'high');
    
    %suggestedfiltorder = minfac*fix(srate/locutoff);
    %fprintf('Suggested filter order  is %d. ',suggestedfiltorder);
else
    disp('ERPLAB_butter1 error: Men at work!')
    return    
end


% Initialize filtered_eeg
smoothdata = zeros(numchan,frames);

hflt = waitbar(0,'Filtering...');
disp('Filtering input data...')
for i = 1:numchan
    smoothdata(i,:) = filtfilt(b,a, double(data(i,:))); % change to filtfilt for non-causal
    waitbar(i/numchan)
end
close(hflt)

