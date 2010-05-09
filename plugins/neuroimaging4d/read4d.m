% read4d() - read 4D Binary datafiles
%	      Return header info, MEG and EEG data
%
% Usage:
%   >> [head, TrialData] = read4d(filename, dataChunks)
%
% Required Input:
%   filename = 4D set data filename
%
% Optional Input:
%   dataChunks = vector containing desired frame numbers(for unsegmented
%                datafiles) or segments (for segmented files). If this
%                input is empty or is not provided then all data will be
%                returned.
% 
% Outputs:
%   head = struct containing header info (see read4dhdr() )
%   MegData = MEG channel data
%
% Author: Christian Wienbruch,
%
% See also: read4dhdr()

function  [head, MegData, EventData] = read4d(filename, epochnumber)

if nargin <1 | nargin >2,
    help read4d;
    return;
end
disp(nargout)
if nargout < 2 | nargout > 3,
	error('2 output args required');
end

% get our header structure
fprintf('Importing binary EGI data file ...\n');
head=read4dhdr(filename);

disp(filename);
datafile = filename(1:length(filename)-4);
[fid,message] = fopen(datafile,'rb','b');
disp(datafile);
if(head.fileformat == 1)
 TrialData=fread(fid, [head.TotalChannels   head.TotalEpochs*head.TotalSlices] ,'short');
end
if(head.fileformat == 2)
 TrialData=fread(fid, [head.TotalChannels head.TotalEpochs*head.TotalSlices] ,'long');
end
if(head.fileformat == 3)
 TrialData=fread(fid, [head.TotalChannels head.TotalEpochs*head.TotalSlices] ,'float');
end
if(head.fileformat == 4)
 disp('data format DOUBLE not yet implemented');
end

MegData = 1.0e15*TrialData(head.MegIndexArray,:);
EegData = 1.0e6*TrialData(head.EegIndexArray,:);
TrigData = TrialData(head.TriggerIndex,:);
RespData = TrialData(head.ResponseIndex,:);
disp(head.TriggerIndex);
disp(head.ResponseIndex);
disp(head.TotalEvents);
for i=1:head.TotalEvents
	EventData(head.Events(i)) = head.EventCodes(i);
end
fclose(fid);
