% eeg_eventhist() - return or plot histogram of event or urevent field values. 
%                   If NO output args, plots the histogram. If the field values 
%                   are not numbers or strings, no histogram is computed.
% Usage:
%         >> figure; eeg_eventhist(EEG.event,'field',bins); % plot histogram
%         >> [fldvals] = eeg_eventhist(EEG.event,'field');  % return field values 
%         >> [fldvals,binNs,binedges] = eeg_eventhist(EEG.event,'field',bins); 
% Inputs:
%
%  Event  - an EEGLAB EEG.event or EEG.urevent structure 
% 'field' - string containing the name of a field in the input Event structure
%  bins   - optional number of bins to use, else vector of bin edges {default: 10}
%           If event field values are strings, this argument is ignored.
% Outputs:
%
% fldvals  - numeric, struct, or cell array vector of field values for each event
%            in the input event order (with [] values, if any, replaced by NaN's or ' 's).
% binNs    - numbers of events in the histogram bins
% binedges - if numeric values, a one-column matrix giving bin edges of the [low,high) bins
%            Else, if string values, cell array containing the string associated with each bin.
% Example:
%         >> [vals,histNs,bins] = eeg_eventhist(EEG.event,'type');
%         %
%         % Returns cell array of event-type strings, numbers of each event type, 
%         % and event type strings, in alphabetic order. No bar() plot produced.
%
% See also:  pop_eventstat(), signalstat(), pop_signalstat().
%
% Author: Scott Makeig, SCCN, Institute for Neural Computation, UCSD, March 26, 2004
%

% 8-20-05 replace found numeric field values [] with NaN to avoid bug -sm
%         replace bin numbers in plot with bin labels if strings; add plot title

function [vals,histNs,outbins] = eeg_eventhist(Event,field,bins)

if nargin < 2
  help eeg_eventhist
  return
end

if nargin < 3
  bins = 10;
end

if isempty(Event)
  error('Event structure is empty');
end
if ~isfield(Event,field)
  error('named field is not an Event field')
end

idx = 0; fld = [];
while isempty(fld) 
   idx = idx+1;
   if idx > length(Event)
     error('All named event fields are empty');
   end
   fld = getfield(Event(idx),field);
   if ischar(fld)
      IS_CHAR = 1;
   elseif isstruct(fld)
      IS_STRUCT = 1;
   elseif ~isempty(fld)
	IS_NUM = 1;
   end
end

if exist('IS_NUM')
  vals = zeros(length(Event),1);
  fprintf('Assuming ''%s'' field values are numeric.\n',field);
elseif exist('IS_CHAR')
  vals = cell(length(Event),1);
  fprintf('Assuming ''%s'' field values are strings.\n',field);
elseif exist('IS_STRUCT')
  vals = repmat(field1,length(Event),1);
  fprintf('Assuming ''%s'' field values are structures.\n',field);
else
  error('Cannot determine field value type')
end

if exist('IS_NUM')
  for k=1:length(Event)
     v = getfield(Event(k),field);
     if isempty(v)
        v = NaN;
     end
     vals(k) = v;
  end
else
  for k=1:length(Event)
     vals{k} = getfield(Event(k),field);
  end
end
if nargout == 1 | exist('IS_STRUCT')
  return      % return vals only, no histogram 
end
%
if exist('IS_NUM')  %%%%%%%%%%%%%%%% numeric values histogram %%%%%%%%%%%%%%%%%%%%%%
%   
  if numel(bins) == 1
    if bins < 3
       error('number of bins must be > 2');
    end

    mn = mean(vals);
    stdev = std(vals);

    binsout = zeros(bins,1);
    fl2 = floor(bins/2);
    for k = -1*fl2:ceil(bins/2)
       binsout(k+fl2+1) = mn+k*stdev;
    end
    binsout(1) = -inf;
    binsout(end) = inf;

    histNs = histc(vals,binsout);
    histNs = histNs(1:end-1);

  else  % accomodate specified bin edges
   histNs = histc(vals,bins);
   histNs = histNs(1:end-1);
  end

  outbins = binsout;
  if nargout == 0
    h = bar(histNs,'histc');
  end
%
else  % exist('IS_CHAR')   %%%%%%%%%%%% string values histogram %%%%%%%%%%%%%%%%%%%%%%
%
   for v=1:length(vals)
     if isempty(cell2mat(vals(v)))
         vals(v) = {' '};
     end
   end
   outbins = unique_bc(vals);
   histNs = zeros(length(outbins)-1,1);
   for k=1:length(outbins)
     histNs(k) = sum(ismember(vals,outbins{k}));
   end

   if nargout == 0
      bar(histNs,1);
      if IS_CHAR
         set(gca,'xticklabel',outbins); % ??? NEEDS MORE WORK - CANT TEST FROM HOME
         yl = get(gca,'ylim');
         set(gca,'ylim',[yl(1) yl(2)*1.1]);
       
         if strcmp(field,'type')
             tl=title(['Histogram of event types']);
         else
             tl=title(['Histogram of event field ''' field ''' values']);
         end
      end
   end
end

return
  
