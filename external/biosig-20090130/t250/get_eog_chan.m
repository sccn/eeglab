function chEOG=get_eog_chan(fn); 
% GET_EOG_CHAN returns EOG channels 
%   
%  EOGchan = GET_EOG_CHAN(...) 

if ischar(fn), 
	HDR=sopen(fn); 
	HDR=sclose(HDR); 
elseif isstruct(fn)
	HDR=fn; 
end; 	

if any(any(isnan(HDR.THRESHOLD(:,1:2))))
	warning('missing threshold(s)'); 
end; 	


%HDR.HIS = histo2(s);
%HDR.RES = hist2res2(HDR.HIS); 
%HDR.RES = y2res(s); 
%HDR.Threshold = repmat([-1,1]/1e2,HDR.NS,1);
HDR.datatype = 'qc:histo'; 

v1 = strmatch('eogv1',lower(HDR.Label));
v2 = strmatch('eogv2',lower(HDR.Label));
v0 = strmatch('eogv',lower(HDR.Label));
v3 = strmatch('eogvp',lower(HDR.Label));
v4 = strmatch('eogvn',lower(HDR.Label));

h1 = strmatch('eogh1',lower(HDR.Label));
h2 = strmatch('eogh2',lower(HDR.Label));
h0 = strmatch('eogh' ,lower(HDR.Label));
h3 = strmatch('eoghp',lower(HDR.Label));
h4 = strmatch('eoghn',lower(HDR.Label));

if length(v1) & length(v2)
	chEOG([v1,v2],1) = [1;-1];
elseif length(v3) & length(v4)
	chEOG([v3,v4],1) = [1;-1];
elseif length(v0)
	chEOG(v0,1) = 1;
end; 

if length(h1) & length(h2)
	chEOG([h1,h2],2) = [1;-1];
elseif length(h3) & length(h4)
	chEOG([h3,h4],2) = [1;-1];
elseif length(h0)
	chEOG(h0,2) = 1;
end; 
[i,j,v]=find(chEOG); 
CHAN = 1:HDR.NS; 
CHAN(i)=[];

if size(chEOG,2)<2, 
	warning('EOG missing')
end; 




