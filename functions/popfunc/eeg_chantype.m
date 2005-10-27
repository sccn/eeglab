function indices = eeg_chantype(data,chantype)
% Returns the channel indices of the desired channel type.

if ischar(chantype), chantype = cellstr(chantype); end
if ~iscell(chantype), error('chantype must be cell array, e.g. {''EEG'', ''EOG''}, or single character string, e.g.''EEG''.'); end

% Define 'datatype' variable, listing the type of each channel.
if isfield(data,'type')
    datatype = {data.type};
elseif isfield(data,'chanlocs') && isfield(data.chanlocs,'type')
    datatype = {data.chanlocs.type};
else error('Incorrect ''data'' input. Should be ''EEG'' or ''loc_file'' structure variable in the format associated with EEGLAB.');
end


k = 1;
for i = 1:length(chantype)
    for j = 1:length(datatype)
        if strcmpi(chantype{i},char(datatype(j)))
            plotchans(k) = j;
            k = k + 1;
        end
    end
end
indices = sort(plotchans);
