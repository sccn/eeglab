function [dat] = fetch_data(data, varargin)

% FETCH_DATA mimics the behaviour of READ_DATA, but for a FieldTrip
% raw data structure instead of a file on disk.
%
% Use as
%   [event] = fetch_data(data, ...)
%
% See also READ_DATA, FETCH_HEADER, FETCH_EVENT

% Copyright (C) 2008, Esther Meeuwissen
%
% $Log: not supported by cvs2svn $
% Revision 1.2  2009/07/06 09:41:18  jansch
% multiple changes. allowing for selection of rpt in frequency data when input
% data has rpttap. allowing for grandaveraging functionality in the case of
% multiple inputs with the same dimensionalities. this is equivalent to the
% XXXgrandaverage functions with keepindividual = 'yes'.
%
% Revision 1.1  2008/11/13 09:55:36  roboos
% moved from fieldtrip/private, fileio or from roboos/misc to new location at fieldtrip/public
%
% Revision 1.2  2008/09/29 21:12:39  roboos
% cleaned up the code from Esther, added copyrights, updated documentation
%

% check whether input is data
data = checkdata(data, 'datatype', 'raw');

% get the options
hdr           = keyval('header',        varargin);
begsample     = keyval('begsample',     varargin);
endsample     = keyval('endsample',     varargin);
chanindx      = keyval('chanindx',      varargin);

if isempty(hdr)
  hdr = fetch_header(data);
end

if isempty(begsample) || isempty(endsample)
  error('begsample and endsample must be specified');
end

if isempty(chanindx)
  chanindx = 1:hdr.nChans;
end

% get trial definition according to original data file
trl    = findcfg(data.cfg, 'trl');
trlnum = length(data.trial);
trllen = zeros(trlnum,1);
for trllop=1:trlnum
  trllen(trllop) = size(data.trial{trllop},2);
end

% check whether data.trial is consistent with trl
if size(trl,1)~=length(data.trial)
  error('trial definition is not internally consistent')
elseif any(trllen~=(trl(:,2)-trl(:,1)+1))
  error('trial definition is not internally consistent')
end

minchan = min(chanindx);
maxchan = max(chanindx);
if minchan<1 || maxchan>hdr.nChans
  error('selected channels are not present in the data')
end

% these are for bookkeeping
maxsample = max(trl(:,2)); 
count     = zeros(1, maxsample);
trialnum  = NaN(1, maxsample);
samplenum = NaN(1, maxsample);

for trllop=1:trlnum
  % make vector with 0= no sample of old data, 1= one sample of old data, 2= two samples of old data, etc
  count(trl(trllop,1):trl(trllop,2)) = count(trl(trllop,1):trl(trllop,2)) + 1;
  % make vector with 1's if samples belong to trial 1, 2's if samples belong to trial 2 etc. overlap/ no data --> Nan
  trialnum(trl(trllop,1):trl(trllop,2)) = trllop;
  % make samplenum vector with samplenrs for each sample in the old trials
  samplenum(trl(trllop,1):trl(trllop,2)) = 1:trllen(trllop);
end

% overlap --> NaN
%trialnum(count>1)  = NaN;
%samplenum(count>1) = NaN;

% make a subselection for the desired samples
count     = count(begsample:endsample);
trialnum  = trialnum(begsample:endsample);
samplenum = samplenum(begsample:endsample);

% check if all samples are present and are not present twice or more 
if any(count==0)
  error('not all requested samples are present in the data');
elseif any(count>1)
%  error('some of the requested samples occur twice in the data');
end

% construct the output data array
dat = zeros(length(chanindx), length(samplenum));
for smplop=1:length(samplenum)
   dat(:, smplop) = data.trial{trialnum(smplop)}(chanindx,samplenum(smplop)); 
end

