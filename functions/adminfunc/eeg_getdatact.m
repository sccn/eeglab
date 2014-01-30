% eeg_getdatact() - get EEG data from a specified dataset or
%                  component activity
%
% Usage:
%       >> signal = eeg_getdatact( EEG );
%       >> signal = eeg_getdatact( EEG, 'key', 'val');
%
% Inputs:
%   EEG       - Input dataset
%
% Optional input:
%   'channel'   - [integer array] read only specific channels.
%                 Default is to read all data channels.
%   'component' - [integer array] read only specific components
%   'projchan'  - [integer or cell array] channel(s) onto which the component
%                 should be projected.
%   'rmcomps'   - [integer array] remove selected components from data
%                 channels. This is only to be used with channel data not
%                 when selecting components.
%   'trialindices' - [integer array] only read specific trials. Default is
%                 to read all trials.
%   'samples'   - [integer array] only read specific samples. Default is
%                 to read all samples.
%   'reshape'   - ['2d'|'3d'] reshape data. Default is '3d' when possible.
%   'verbose'   - ['on'|'off'] verbose mode. Default is 'on'.
%
% Outputs:
%   signal      - EEG data or component activity
%
% Author: Arnaud Delorme, SCCN & CERCO, CNRS, 2008-
%
% See also: eeg_checkset()

% Copyright (C) 15 Feb 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [data boundaries] = eeg_getdatact( EEG, varargin);

data = [];
if nargin < 1
    help eeg_getdatact;
    return;
end;

% reading data from several datasets and concatening it
% -----------------------------------------------------
if iscell(EEG) || (~isstr(EEG) && length(EEG) > 1)
    % decode some arguments
    % ---------------------
    trials  = cell(1,length(EEG));
    rmcomps = cell(1,length(EEG));
    for iArg = length(varargin)-1:-2:1
        if strcmpi(varargin{iArg}, 'trialindices')
            trials = varargin{iArg+1};
            varargin(iArg:iArg+1) = [];
        elseif strcmpi(varargin{iArg}, 'rmcomps')
            rmcomps = varargin{iArg+1};
            varargin(iArg:iArg+1) = [];
        end;
    end;
    if isnumeric(rmcomps), rmtmp = rmcomps; rmcomps = cell(1,length(EEG)); rmcomps(:) = { rmtmp }; end;
        
    % concatenate datasets
    % --------------------
    data = [];
    boundaries = [];
    for dat = 1:length(EEG)
        if iscell(EEG)
             [tmpdata datboundaries] = eeg_getdatact(EEG{dat}, 'trialindices', trials{dat}, 'rmcomps', rmcomps{dat}, varargin{:} );
        else [tmpdata datboundaries] = eeg_getdatact(EEG(dat), 'trialindices', trials{dat}, 'rmcomps', rmcomps{dat}, varargin{:} );
        end;
        if isempty(data), 
            data       = tmpdata;
            boundaries = datboundaries;
        else
            if all([ EEG.trials ] == 1) % continuous data
                if size(data,1) ~= size(tmpdata,1), error('Datasets to be concatenated do not have the same number of channels'); end;

                % adding boundaries
                if ~isempty(datboundaries)
                    boundaries = [boundaries datboundaries size(data,2)];
                else
                    boundaries = [boundaries size(data,2)];
                end;
                data(:,end+1:end+size(tmpdata,2)) = tmpdata; % concatenating trials
            else
                if size(data,1) ~= size(tmpdata,1), error('Datasets to be concatenated do not have the same number of channels'); end;
                if size(data,2) ~= size(tmpdata,2), error('Datasets to be concatenated do not have the same number of time points'); end;
                data(:,:,end+1:end+size(tmpdata,3)) = tmpdata; % concatenating trials
            end;
        end;
    end;
    return;
end;

% if string load dataset
% ----------------------
if isstr(EEG)
    EEG = pop_loadset('filename', EEG, 'loadmode', 'info');
end;

opt = finputcheck(varargin, { ...
    'channel'   'integer' {} [];
    'verbose'   'string'  { 'on','off' } 'on';
    'reshape'   'string'  { '2d','3d' }  '3d';
    'projchan'  {'integer','cell' } { {} {} } [];
    'component' 'integer' {} [];
    'samples'   'integer' {} [];
    'interp'    'struct'  { }        struct([]);
    'trialindices' {'integer','cell'} { {} {} } [];
    'rmcomps'      {'integer','cell'} { {} {} } [] }, 'eeg_getdatact');

if isstr(opt), error(opt); end;
channelNotDefined = 0;
if isempty(opt.channel), opt.channel = [1:EEG.nbchan]; channelNotDefined = 1;
elseif isequal(opt.channel, [1:EEG.nbchan]) && ~isempty(opt.interp) channelNotDefined = 1;
end;
if isempty(opt.trialindices), opt.trialindices = [1:EEG.trials]; end;
if iscell( opt.trialindices), opt.trialindices = opt.trialindices{1}; end;
if iscell(opt.rmcomps     ), opt.rmcomps      = opt.rmcomps{1};      end;
if (~isempty(opt.rmcomps) | ~isempty(opt.component)) & isempty(EEG.icaweights)
    error('No ICA weight in dataset');
end;

if strcmpi(EEG.data, 'in set file')
    EEG = pop_loadset('filename', EEG.filename, 'filepath', EEG.filepath);
end;

% get data boundaries if continuous data
% --------------------------------------
boundaries = [];
if nargout > 1 && EEG.trials == 1 && ~isempty(EEG.event) && isfield(EEG.event, 'type') && isstr(EEG.event(1).type)
    if ~isempty(opt.samples)
        disp('WARNING: eeg_getdatact.m, boundaries are not accurate when selecting data samples');
    end;
    tmpevent = EEG.event;
    tmpbound = strmatch('boundary', lower({ tmpevent.type }));
    if ~isempty(tmpbound)
        boundaries = [tmpevent(tmpbound).latency ]-0.5;
    end;
end;

% getting channel or component activation
% ---------------------------------------
filename = fullfile(EEG.filepath, [ EEG.filename(1:end-4) '.icaact' ] );
if ~isempty(opt.component) & ~isempty(EEG.icaact)
    
    data = EEG.icaact(opt.component,:,:);
    
elseif ~isempty(opt.component) & exist(filename)
    
    % reading ICA file
    % ----------------
    data = repmat(single(0), [ length(opt.component) EEG.pnts EEG.trials ]);
    fid = fopen( filename, 'r', 'ieee-le'); %little endian (see also pop_saveset)
    if fid == -1, error( ['file ' filename ' could not be open' ]); end;
    for ind = 1:length(opt.component)
        fseek(fid, (opt.component(ind)-1)*EEG.pnts*EEG.trials*4, -1);
        data(ind,:) = fread(fid, [EEG.trials*EEG.pnts 1], 'float32')';
    end;
    fclose(fid);
    
elseif ~isempty(opt.component)
    
    if isempty(EEG.icaact)
        data = eeg_getdatact( EEG );
        data = (EEG.icaweights(opt.component,:)*EEG.icasphere)*data(EEG.icachansind,:);
    else
        data = EEG.icaact(opt.component,:,:);
    end;
    
else
    if isnumeric(EEG.data) % channel
    
        data = EEG.data;
    
    else % channel but no data loaded

        filename = fullfile(EEG.filepath, EEG.data);
        fid = fopen( filename, 'r', 'ieee-le'); %little endian (see also pop_saveset)
        if fid == -1
            error( ['file ' filename ' not found. If you have renamed/moved' 10 ...
                'the .set file, you must also rename/move the associated data file.' ]);
        else
            if strcmpi(opt.verbose, 'on')
                fprintf('Reading float file ''%s''...\n', filename);
            end;
        end;

        % old format = .fdt; new format = .dat (transposed)
        % -------------------------------------------------
        datformat = 0;
        if length(filename) > 3
            if strcmpi(filename(end-2:end), 'dat')
                datformat = 1;
            end;
        end;
        EEG.datfile = EEG.data;

        % reading data file
        % -----------------
        eeglab_options;
        if length(opt.channel) == EEG.nbchan && option_memmapdata
            fclose(fid);
            data = mmo(filename, [EEG.nbchan EEG.pnts EEG.trials], false);
            %data = memmapdata(filename, [EEG.nbchan EEG.pnts EEG.trials]);
        else
            if datformat
                if length(opt.channel) == EEG.nbchan || ~isempty(opt.interp)
                    data = fread(fid, [EEG.trials*EEG.pnts EEG.nbchan], 'float32')';
                else
                    data = repmat(single(0), [ length(opt.channel) EEG.pnts EEG.trials ]);
                    for ind = 1:length(opt.channel)
                        fseek(fid, (opt.channel(ind)-1)*EEG.pnts*EEG.trials*4, -1);
                        data(ind,:) = fread(fid, [EEG.trials*EEG.pnts 1], 'float32')';
                    end;
                    opt.channel = [1:size(data,1)];
                end;
            else
                data = fread(fid, [EEG.nbchan Inf], 'float32');
            end;
            fclose(fid);
        end;
        
    end;
    
    % subracting components from data
    % -------------------------------
    if ~isempty(opt.rmcomps)
        if strcmpi(opt.verbose, 'on')
            fprintf('Removing %d artifactual components\n', length(opt.rmcomps));
        end;
        rmcomps = eeg_getdatact( EEG, 'component', opt.rmcomps); % loaded from file
        rmchan    = [];
        rmchanica = [];
        for index = 1:length(opt.channel)
            tmpicaind = find(opt.channel(index) == EEG.icachansind);
            if ~isempty(tmpicaind)
                rmchan    = [ rmchan    index ];
                rmchanica = [ rmchanica tmpicaind ];
            end;
        end;
        data(rmchan,:) = data(rmchan,:) - EEG.icawinv(rmchanica,opt.rmcomps)*rmcomps(:,:);
        
        %EEG = eeg_checkset(EEG, 'loaddata');
        %EEG = pop_subcomp(EEG, opt.rmcomps);
        %data = EEG.data(opt.channel,:,:);
        
        %if strcmpi(EEG.subject, 'julien') & strcmpi(EEG.condition, 'oddball') & strcmpi(EEG.group, 'after')
        %    jjjjf
        %end;
    end;
    
    if ~isempty(opt.interp)
        EEG.data   = data;
        EEG.event  = [];
        EEG.epoch  = [];
        EEG = eeg_interp(EEG, opt.interp, 'spherical');
        data = EEG.data;
        if channelNotDefined, opt.channel = [1:EEG.nbchan]; end;
    end;

    if ~isequal(opt.channel, [1:EEG.nbchan])
        data = data(intersect(opt.channel,[1:EEG.nbchan]),:,:);
    end;
end;


% projecting components on data channels
% --------------------------------------
if ~isempty(opt.projchan)
    if iscell(opt.projchan)
        opt.projchan = std_chaninds(EEG, opt.projchan);
    end;
    
    finalChanInds = [];
    for iChan = 1:length(opt.projchan)
        tmpInd = find(EEG.icachansind == opt.projchan(iChan));
        if isempty(tmpInd)
            error(sprintf('Warning: can not backproject component on channel %d (not used for ICA)\n', opt.projchan(iChan)));
        end;
        finalChanInds = [ finalChanInds tmpInd ];
    end;
    
    data = EEG.icawinv(finalChanInds, opt.component)*data(:,:);
end;

if size(data,2)*size(data,3) ~= EEG.pnts*EEG.trials
    disp('WARNING: The file size on disk does not correspond to the dataset, file has been truncated');
end;
try,
    if EEG.trials == 1, EEG.pnts = size(data,2); end;
    if  strcmpi(opt.reshape, '3d')
         data = reshape(data, size(data,1), EEG.pnts, EEG.trials);
    else data = reshape(data, size(data,1), EEG.pnts*EEG.trials);
    end;
catch
    error('The file size on disk does not correspond to the dataset information.');
end;

% select trials
% -------------
if length(opt.trialindices) ~= EEG.trials
    data = data(:,:,opt.trialindices);
end;
if ~isempty(opt.samples)
    data = data(:,opt.samples,:);
end;
