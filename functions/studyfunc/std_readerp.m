% std_readerp() - load ERP measures for data channels or 
%                  for all components of a specified cluster.
%                  Called by plotting functions
%                  std_envtopo(), std_erpplot(), std_erspplot(), ...
% Usage:
%         >> [STUDY, datavals, times, setinds, cinds] = ...
%                   std_readerp(STUDY, ALLEEG, varargin);
% Inputs:
%       STUDY - studyset structure containing some or all files in ALLEEG
%      ALLEEG - vector of loaded EEG datasets
%
% Optional inputs:
%  'design'    - [integer] read files from a specific STUDY design. Default
%                is empty (use current design in STUDY.currentdesign).
%  'channels'  - [cell] list of channels to import {default: none}
%  'clusters'  - [integer] list of clusters to import {[]|default: all but
%                the parent cluster (1) and any 'NotClust' clusters}
%  'singletrials' - ['on'|'off'] load single trials spectral data (if 
%                available). Default is 'off'.
%  'subject'   - [string] select a specific subject {default:all}
%  'component' - [integer] select a specific component in a cluster
%                 {default:all}
%
% ERP specific optional inputs:
%  'timerange' - [min max] time range {default: whole measure range}
%  'componentpol' - ['on'|'off'] invert ERP component sign based on
%                   scalp map match with component scalp map centroid.
%                   {default:'on'}
%
% Output:
%  STUDY    - updated studyset structure
%  datavals  - [cell array] erp data (the cell array size is 
%             condition x groups)
%  times    - [float array] array of time values
%  setinds  - [cell array] datasets indices
%  cinds    - [cell array] channel or component indices
%
% Example:
%  std_precomp(STUDY, ALLEEG, { ALLEEG(1).chanlocs.labels }, 'erp', 'on');
%  [erp times] = std_readerp(STUDY, ALLEEG, 'channels', { ALLEEG(1).chanlocs(1).labels });
%
% Author: Arnaud Delorme, CERCO, 2006-

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function [STUDY, datavals, xvals, setinds, allinds] = std_readerp(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readerp;
    return;
end
if ~isstruct(ALLEEG) 
    % old calling format
    % ------------------
    EEG = STUDY(ALLEEG);
    filename   = fullfile(EEG.filepath, EEG.filename(1:end-4));
    comporchan = varargin{1};
    options = {'measure', 'erp'};
    if length(varargin) > 1, options = { options{:} 'timelimits', varargin{2} }; end;
    if comporchan(1) > 0
        [datavals tmp xvals] = std_readfile(filename, 'components',comporchan, options{:});
    else
       [datavals tmp xvals] = std_readfile(filename, 'channels', comporchan, options{:});
    end;
    STUDY    = datavals';
    datavals = xvals;
    return;
end;

STUDY = pop_erpparams(STUDY, 'default');
STUDY = pop_specparams(STUDY, 'default');
[opt moreopts] = finputcheck( varargin, { ...
    'design'        'integer' []             STUDY.currentdesign;
    'channels'      'cell'    []             {};
    'clusters'      'integer' []             [];
    'timerange'     'real'    []             STUDY.etc.erpparams.timerange;
    'freqrange'     'real'    []             STUDY.etc.specparams.freqrange;
    'datatype'      'string'  { 'erp','spec' } 'erp';
    'rmsubjmean'    'string'  { 'on','off' } 'off';
    'singletrials'  'string'  { 'on','off' } 'off';
    'componentpol'  'string'  { 'on','off' } 'on';
    'component'     'integer' []             [];
    'subject'       'string'  []             '' }, ...
    'std_readerp', 'ignore');
if isstr(opt), error(opt); end;
nc = max(length(STUDY.design(opt.design).variable(1).value),1);
ng = max(length(STUDY.design(opt.design).variable(2).value),1);
paired1 = STUDY.design(opt.design).variable(1).pairing;
paired2 = STUDY.design(opt.design).variable(2).pairing;

dtype = opt.datatype;

% find channel indices
% --------------------
if ~isempty(opt.channels)
     allChangrp = lower({ STUDY.changrp.name });
     finalinds = std_chaninds(STUDY, opt.channels);
else finalinds = opt.clusters;
end;

for ind = 1:length(finalinds) % scan channels or components

    % find indices
    % ------------
    if ~isempty(opt.channels)
        tmpstruct = STUDY.changrp(finalinds(ind));
        allinds       = tmpstruct.allinds;
        setinds       = tmpstruct.setinds;
    else
        tmpstruct = STUDY.cluster(finalinds(ind));
        allinds       = tmpstruct.allinds;
        setinds       = tmpstruct.setinds;
    end;

    % check if data is already here
    % -----------------------------
    dataread = 0;
    if strcmpi(dtype, 'erp'), eqtf = isequal( STUDY.etc.erpparams.timerange, opt.timerange);
    else                      eqtf = isequal(STUDY.etc.specparams.freqrange, opt.freqrange) && ...
                                     isequal( STUDY.etc.specparams.subtractsubjectmean, opt.rmsubjmean);
    end;
    if strcmpi(opt.singletrials,'off')
        if isfield(tmpstruct, [ dtype 'data' ]) && ~isempty(getfield(tmpstruct, [ dtype 'data' ])) && eqtf
            dataread = 1;
        end;
    else
        if isfield(tmpstruct, [ dtype 'datatrials' ]) && eqtf
            tmpdat = getfield(tmpstruct, [ dtype 'datatrials' ]);
            range  = fastif( strcmpi(dtype, 'erp'), 'erptimes', 'specfreqs');
            if isfield(tmpstruct, range)
                if ~isempty(opt.channels) && ~isempty(tmpdat) && strcmpi(getfield(tmpstruct, [ dtype 'trialinfo' ]), opt.subject)
                    if size(tmpdat{1},2) == length(getfield(tmpstruct, range))
                        dataread = 1; 
                    end;
                elseif isempty(opt.channels) && isequal(getfield(tmpstruct, [ dtype 'trialinfo' ]), opt.component) 
                    if size(tmpdat{1},2) == length(getfield(tmpstruct, range))
                        dataread = 1; 
                    end;
                end;
            end;
        end;
    end;
    
    if ~dataread
        % reserve arrays
        % --------------
        alldata        = cell( nc, ng );
        setinfoIndices = cell( nc, ng );
        tmpind  = 1; while(isempty(setinds{tmpind})), tmpind = tmpind+1; end;
        setinfo = STUDY.design(opt.design).cell;
        tmpchanlocs = ALLEEG(setinfo(1).dataset(1)).chanlocs;
        chanlab = { tmpchanlocs.labels };
        nonemptyindex = ~cellfun(@isempty, allinds);
        nonemptyindex = find(nonemptyindex(:));
        optGetparams = { 'measure', dtype, 'getparamonly', 'on', 'singletrials', opt.singletrials, 'timelimits', opt.timerange, 'freqlimits', opt.freqrange, 'setinfoinds', 1};
        if ~isempty(opt.channels), [ tmp params xvals] = std_readfile(setinfo(setinds{nonemptyindex(1)}(1)), optGetparams{:}, 'channels'  , allChangrp(allinds{nonemptyindex(1)}(1)));
        else                       [ tmp params xvals] = std_readfile(setinfo(setinds{nonemptyindex(1)}(1)), optGetparams{:}, 'components', allinds{nonemptyindex(1)}(1));
        end;
        
        % read the data and select channels
        % ---------------------------------
        fprintf([ 'Reading ' dtype ' data...' ]);
        if strcmpi(dtype, 'erp'), opts = { 'timelimits', opt.timerange };
        else                      opts = { 'freqlimits', opt.freqrange };
        end;
        if strcmpi(opt.singletrials, 'on')
            if strcmpi(params.singletrials, 'off')
                fprintf('\n');
                errordlg2('No single trial data - recompute data files');
                datavals = [];
                return;
            end;            
            opts        = { opts{:} 'singletrials' 'on' };
            for c = 1:nc
                for g = 1:ng
                    if ~isempty(setinds{c,g})
                        if ~isempty(opt.channels), [alldata{c, g} z z z z setinfoIndices{c, g}] = std_readfile( setinfo(setinds{c,g}(:)), 'measure', dtype, opts{:}, 'channels'  , opt.channels(ind), 'setinfoinds', setinds{c,g}(:));
                        else                       [alldata{c, g} z z z z setinfoIndices{c, g}] = std_readfile( setinfo(setinds{c,g}(:)), 'measure', dtype, opts{:}, 'components', allinds{c,g},      'setinfoinds', setinds{c,g}(:));
                        end;
                    end;
                end;
            end;
        else
            for c = 1:nc
                for g = 1:ng
                    if ~isempty(setinds{c,g})
                        if ~isempty(opt.channels), [alldata{c, g}] = std_readfile( setinfo(setinds{c,g}(:)), 'measure', dtype, opts{:}, 'channels'  , opt.channels(ind));
                        else                       [alldata{c, g}] = std_readfile( setinfo(setinds{c,g}(:)), 'measure', dtype, opts{:}, 'components', allinds{c,g});
                        end;
                    end;
                end;
            end;
        end;
        fprintf('\n');
        
        % inverting component polaritites
        % -------------------------------
        if isempty(opt.channels) && strcmpi(dtype, 'erp')
            if strcmpi(opt.singletrials, 'on')
                disp('Warning: component ERP polarity cannot (yet) be inverted for single trials');
            else
                STUDY = std_readtopoclust(STUDY, ALLEEG, finalinds(ind));
                if isfield(STUDY.cluster, 'topopol') && ~isempty(STUDY.cluster(finalinds(ind)).topopol)
                    [ tmpstruct tmp1 tmp2 topopolcell] = std_setcomps2cell(STUDY, STUDY.cluster(finalinds(ind)).sets, STUDY.cluster(finalinds(ind)).comps, STUDY.cluster(finalinds(ind)).topopol);
                    disp('Inverting ERP component polarities based on scalp map polarities');
                    for index = 1:length(alldata(:))
                        for comps = 1:size(alldata{index},2)
                            alldata{index}(:,comps) = alldata{index}(:,comps)*topopolcell{index}(comps);
                        end;
                    end;
                else
                    disp('Cluster topographies absent - cannot adjust single component ERP polarities');
                end;
            end;
        end;

        % remove mean of each subject across groups and conditions
        if strcmpi(dtype, 'spec') && strcmpi(opt.rmsubjmean, 'on') && ~isempty(opt.channels)
            disp('Removing subject''s average spectrum based on pairing settings');
            if strcmpi(paired1, 'on') && strcmpi(paired2, 'on') && (nc > 1 || ng > 1)
                disp('Removing average spectrum for both indep. variables');
                meanpowbase = computemeanspectrum(alldata(:), opt.singletrials);
                alldata     = removemeanspectrum(alldata, meanpowbase);
            elseif strcmpi(paired1, 'on') && ng > 1
                disp('Removing average spectrum for first indep. variables (second indep. var. is unpaired)');
                for g = 1:ng        % ng = number of groups
                    meanpowbase  = computemeanspectrum(alldata(:,g), opt.singletrials);
                    alldata(:,g) = removemeanspectrum(alldata(:,g), meanpowbase);
                end;
            elseif strcmpi(paired2, 'on') && nc > 1
                disp('Removing average spectrum for second indep. variables (first indep. var. is unpaired)');
                for c = 1:nc        % ng = number of groups
                    meanpowbase  = computemeanspectrum(alldata(c,:), opt.singletrials);
                    alldata(c,:) = removemeanspectrum(alldata(c,:), meanpowbase);
                end;
            else
                disp('Not removing average spectrum baseline (both indep. variables are unpaired');
            end;
        end;
        
        if strcmpi(opt.singletrials, 'on')
             tmpstruct = setfield( tmpstruct, [ dtype 'datatrials' ], alldata);
             tmpstruct = setfield( tmpstruct, [ 'setindstrials' ], setinfoIndices);             
             if ~isempty(opt.channels)
                  tmpstruct = setfield( tmpstruct, [ dtype 'trialinfo' ], opt.subject);
             else tmpstruct = setfield( tmpstruct, [ dtype 'trialinfo' ], opt.component);
             end;
        else tmpstruct = setfield( tmpstruct, [ dtype 'data' ], alldata);
        end;
        if strcmpi(dtype, 'spec'), tmpstruct.specfreqs = xvals;
        else                       tmpstruct.erptimes  = xvals;
        end;
        
        % copy results to structure
        % -------------------------
        fieldnames = { [ dtype 'data' ]  [ dtype 'freqs' ] [ dtype 'datatrials' ] ...
                       [ dtype 'times' ] [ dtype 'trialinfo' ] 'allinds' 'setinds' 'setindstrials' };
        for f = 1:length(fieldnames)
            if isfield(tmpstruct, fieldnames{f}),
                tmpdata = getfield(tmpstruct, fieldnames{f});
                if ~isempty(opt.channels)
                     STUDY.changrp = setfield(STUDY.changrp, { finalinds(ind) }, fieldnames{f}, tmpdata);
                else STUDY.cluster = setfield(STUDY.cluster, { finalinds(ind) }, fieldnames{f}, tmpdata);
                end;
            end;
        end;
    end;
end;

% if several channels, agregate them
% ----------------------------------
allinds   = finalinds;
if ~isempty(opt.channels)
    structdat = STUDY.changrp;
    datavals = cell(nc, ng);
    if strcmpi( dtype, 'spec'), xvals = getfield(structdat(allinds(1)), [ dtype 'freqs' ]);
    else                        xvals = getfield(structdat(allinds(1)), [ dtype 'times' ]);
    end;
    for ind =  1:length(datavals(:))
        if strcmpi(opt.singletrials, 'on')
             tmpdat = getfield(structdat(allinds(1)), [ dtype 'datatrials' ]);
        else tmpdat = getfield(structdat(allinds(1)), [ dtype 'data' ]);
        end;
        if ~isempty(tmpdat{ind})
            datavals{ind} = zeros([ size(tmpdat{ind}) length(allinds)]);
            for chan = 1:length(allinds)
                if strcmpi(opt.singletrials, 'on')
                     tmpdat = getfield(structdat(allinds(chan)), [ dtype 'datatrials' ]);
                else tmpdat = getfield(structdat(allinds(chan)), [ dtype 'data' ]);
                end;
                datavals{ind}(:,:,chan) = tmpdat{ind}; % only works for interpolated data
            end;
        
            datavals{ind} = squeeze(permute(datavals{ind}, [1 3 2])); % time elec subjects
        end;
    end;
    setinds  = structdat(allinds(1)).setinds;
    if ~isempty(opt.subject)
        if strcmpi(opt.singletrials, 'on')
            datavals = std_selsubject(datavals, opt.subject, structdat(allinds(1)).setindstrials, { STUDY.design(opt.design).cell.case }, 2); 
        else
            datavals = std_selsubject(datavals, opt.subject, setinds, { STUDY.design(opt.design).cell.case }, 2); 
        end;
    end;
else
    if strcmpi(opt.singletrials, 'on')
         datavals = getfield(STUDY.cluster(allinds(1)), [ dtype 'datatrials' ]);
    else datavals = getfield(STUDY.cluster(allinds(1)), [ dtype 'data' ]);
    end;
    if strcmpi( dtype, 'spec'), xvals = getfield(STUDY.cluster(allinds(1)), [ dtype 'freqs' ]);
    else                        xvals = getfield(STUDY.cluster(allinds(1)), [ dtype 'times' ]);
    end;
    compinds = STUDY.cluster(allinds(1)).allinds;
    setinds  = STUDY.cluster(allinds(1)).setinds;
    if ~isempty(opt.component) && length(allinds) == 1 && strcmpi(opt.singletrials,'off')
        datavals = std_selcomp(STUDY, datavals, allinds, setinds, compinds, opt.component);
    end;
end;

% compute mean spectrum 
% ---------------------
function meanpowbase = computemeanspectrum(spectrum, singletrials)

    try
        len = length(spectrum(:));
        count = 0;
        for index = 1:len
            if ~isempty(spectrum{index})
                if strcmpi(singletrials, 'on')
                    if count == 0, meanpowbase = mean(spectrum{index},2);
                    else           meanpowbase = meanpowbase + mean(spectrum{index},2);
                    end;
                else
                    if count == 0, meanpowbase = spectrum{index};
                    else           meanpowbase = meanpowbase + spectrum{index};
                    end;
                end;
                count = count+1;
            end;
        end;
        meanpowbase = meanpowbase/count;
    catch,
        error([ 'Problem while subtracting mean spectrum.' 10 ...
                'Common spectrum subtraction is performed based on' 10 ...
                'pairing settings in your design. Most likelly, one' 10 ...
                'independent variable should not have its data paired.' ]);
    end;
        
% remove mean spectrum 
% --------------------
function spectrum = removemeanspectrum(spectrum, meanpowbase)
    for g = 1:size(spectrum,2)        % ng = number of groups
        for c = 1:size(spectrum,1)
            if ~isempty(spectrum{c,g}) && ~isempty(spectrum{c,g})
                if size(spectrum{c,g},2) ~= size(meanpowbase, 2)
                     tmpmeanpowbase = repmat(meanpowbase, [1 size(spectrum{c,g},2)]);
                else tmpmeanpowbase = meanpowbase;
                end;
                spectrum{c,g} = spectrum{c,g} - tmpmeanpowbase;
            end;
        end;
    end;

