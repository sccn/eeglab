function  [TrialData] = EGIread(filename)
% [TrialData] = EGIread(filename)
% reads in version 2 and 3 of EGI Simple Binary data files,
% filename = EGI data filename
% TrialData = EEG channel data 
% EventCodes = event codes 
% Samp_Rate = sampling rate
% NChan = #of channels
% NSamp = sampling rate
% Segments = # of epochs
% NEvent = number of events
%

[fid,message] = fopen(filename,'rb','b');
if (fid == 0),
   error(message);
end

version = fread(fid,1,'integer*4');

if (version == 2 | version == 3),
	disp('Reading EGI Simple Binary format data...');
else
	error('EGI Simple Binary Version 2 and 3 supported only.');
end;

year = fread(fid,1,'integer*2');
month = fread(fid,1,'integer*2');
day = fread(fid,1,'integer*2');
hour = fread(fid,1,'integer*2');
minute = fread(fid,1,'integer*2');
second = fread(fid,1,'integer*2');
millisecond = fread(fid,1,'integer*4');
Samp_Rate = fread(fid,1,'integer*2');
NChan = fread(fid,1,'integer*2');
Gain = fread(fid,1,'integer*2');
Bits = fread(fid,1,'integer*2');
Range = fread(fid,1,'integer*2');

Segments = 1;

if (version == 3),
	Categories = fread(fid,1,'integer*2');
	if (Categories),
		for i=1:Categories,
			catname_len(i) = fread(fid,1,'uchar');
			Catname{i} = char(fread(fid,catname_len(i),'uchar'));
		end
	end
	Segments = fread(fid,1,'integer*2');
end 

NSamp = fread(fid,1,'integer*4');
NEvent = fread(fid,1,'integer*2');

if (NEvent == 0)
	EventCodes(1,1:4) = 'none';
else
	for i = 1:NEvent
		EventCodes(i,1:4) = fread(fid,[1,4],'uchar');
	end
end;

FrameVals = NChan + NEvent;
TrialData = zeros(FrameVals,NSamp*Segments);

for i=1:Segments,
	if (version == 3),
		SegmentCatIndex(i) = fread(fid,1,'integer*2');
		SegmentStartTime(i) = fread(fid,1,'integer*4');
	end

	TrialData(:,[1+(i-1)*NSamp:i*NSamp]) = ...
	   fread(fid,[FrameVals,NSamp],'integer*2');
end

fclose(fid);

if (NEvent > 0),
	EventData = TrialData(NChan+1:end,:);
	TrialData = TrialData(1:NChan,:);
end 
