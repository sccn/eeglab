function [data,srate,numSamples,labels,events] = eeg_load_scan4_cnt_data(fid,varargin)

% eeg_load_scan4_cnt_data - Read Neuroscan 4.x Continuous EEG Data
% 
% [data,srate,numSamples,labels,events] = eeg_load_scan4_cnt_data(fid,channels,range)
% 
%   ------------- IN ------------------------------------------------------------                      
%   fid                     ->  file identifier
%   channels (optional)     ->  vector of channel numbers to load or string 'all'.
%                               [default = 'all']                   
%   range (optional)        ->  2-element vector containing start and stop points
%                               (inclusive) or string 'all' [default = 'all']
%
%   ------------- OUT ------------------------------------------------------------
%   data                    <-  matrix of electrodes in columns 
%                               (ascending order, no repeats)
%   srate (optional)        <-  sample rate of data
%   numSamples (optional)   <-  total samples in CNT file
%   labels (optional)       <-  cell array of channel labels
%   events (optional)       <-  matrix of event info:
%                                   column 1 = event stim type
%                                   column 2 = offset in points
%
%   Note: Works only with Scan 4.1+ data files
%
%   See also: eeg_load_scan4_cnt_event
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:52 $

% Licence:  GNU GPL, no express or implied warranties
% History:  2002, Sean.Fitzgibbon@flinders.edu.au
%           06/2002, Darren.Weber_at_radiology.ucsf.edu
%                    adapted to eeg_toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --------------- Load Parameters -----------------------------
fseek(fid,370,'bof');
numChan = fread(fid,1,'ushort');

fseek(fid,376,'bof');
srate = fread(fid,1,'ushort');

fseek(fid,886,'bof');
eventPos = fread(fid,1,'long');

dataPos = 900+(75*numChan);
numSamples = ((eventPos - dataPos)/numChan)/2;

fseek(fid,900,'bof');
labels = cellstr(deblank(char(fread(fid,[10,numChan],'10*char',75-10)')));

fseek(fid,947,'bof');
baseline = fread(fid,numChan,'short',75-2)';

fseek(fid,959,'bof');
sensitivity = fread(fid,numChan,'float',75-4)';

fseek(fid,971,'bof');
calibration = fread(fid,numChan,'float',75-4)';

% --------------- Set Channels & Ranges -----------------

if (nargin == 1)
    channels = 'all';
    range = 'all';
elseif (nargin == 2)
    channels = varargin{1};
    range = 'all';
elseif  (nargin == 3)
    channels = varargin{1};
    range = varargin{2};
end

if isa(range,'char') & (range == 'all')
    start = 0;
    numPoints = numSamples;
    range = [1 numSamples];
else
    start = (range(1)-1)*numChan*2;
    numPoints = range(2)-range(1)+1;
end





% Must check than range is within numSamples!






if isa(channels,'char') & (channels == 'all')
    channels = [1:numChan];
else
    channels = sort(channels);
    index = [];
    for i = 1:length(channels)-1
        if (channels(i) == channels(i+1))
            index = [index i+1];
        end
    end
    channels(index) = [];
end


% --------------- Read Events ---------------------------


fseek(fid,eventPos+1,'bof');
eventSize   = fread(fid,1,'long');
eventOffset = fread(fid,1,'long');

fseek(fid,eventPos + 9 + eventOffset,'bof');
stimType = fread(fid,eventSize/19,'short',19-2);

fseek(fid,eventPos + 9 + eventOffset + 4,'bof');
stimOffset = fread(fid,eventSize/19,'long',19-4);

stimOffset = stimOffset - (900 + (75 * numChan));
stimOffset = stimOffset ./ (2 * numChan);

events = [stimType stimOffset];


% --------------- Read Data -----------------------------

if (length(channels) < numChan)
    data = zeros(numPoints,length(channels));
    for i = 1:length(channels)
        fseek(fid,dataPos+start,'bof');
        fseek(fid,(channels(i)-1)*2,'cof');
        data(:,i) = fread(fid,numPoints,'short',(numChan-1)*2);
    end
elseif (length(channels) == numChan)
    fseek(fid,dataPos+start,'bof');
    data = zeros(numChan,numPoints);
    data = fread(fid,[numChan,numPoints],'short');
    data = data';
end


% --------------- scale to microvolts -----------------------------

scaleFactor = ((sensitivity .* calibration)./204.8);
scaleFactor = scaleFactor(ones(size(data,1),1),channels);

baseline = baseline(ones(size(data,1),1),channels);

data = (data - baseline) .* scaleFactor;

frewind(fid);

return
