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
%                is empty (use current design in STUDY.currentdesign). Use
%                NaN to create a design with with all the data.
%  'channels'  - [cell] list of channels to import {default: none}
%  'clusters'  - [integer] list of clusters to import {[]|default: all but
%                the parent cluster (1) and any 'NotClust' clusters}
%  'singletrials' - ['on'|'off'] load single trials spectral data (if 
%                available). Default is 'off'.
%  'subject'   - [string] select a specific subject {default:all}
%  'component' - [integer] select a specific component in a cluster.
%                This is the index of the component in the cluster not the
%                component number {default:all}
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

function [STUDY, datavals, xvals, yvals, events, params] = std_readerp(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readerp;
    return;
end

STUDY = pop_erpparams(STUDY, 'default');
STUDY = pop_specparams(STUDY, 'default');
STUDY = pop_erspparams(STUDY, 'default');
[opt moreopts] = finputcheck( varargin, { ...
    'design'        'integer' []             STUDY.currentdesign;
    'channels'      'cell'    []             {};
    'clusters'      'integer' []             [];
    'timerange'     'real'    []             [];
    'freqrange'     'real'    []             [];
    'datatype'      'string'  { 'erp','spec' 'ersp' 'itc' } 'erp';
    'rmsubjmean'    'string'  { 'on','off' } 'off';
    'subbaseline'   'string'  { 'on','off' } STUDY.etc.erspparams.subbaseline; % subtract common baseline (ERSP only)
    'singletrials'  'string'  { 'on','off' } 'off';
    'componentpol'  'string'  { 'on','off' } 'on';
    'component'     'integer' []             [];
    'subject'       'string'  []             '' }, ...
    'std_readerp', 'ignore');
if isstr(opt), error(opt); end;

dtype = opt.datatype;

% get the file extension
% ----------------------
tmpDataType = opt.datatype;
if strcmpi(opt.datatype, 'ersp') || strcmpi(opt.datatype, 'itc'), 
    tmpDataType = 'timef'; 
    if isempty(opt.timerange), opt.timerange = STUDY.etc.erpparams.timerange;  end;
    if isempty(opt.freqrange), opt.freqrange = STUDY.etc.specparams.freqrange; end;
else
    if isempty(opt.timerange), opt.timerange = STUDY.etc.erspparams.timerange;  end;
    if isempty(opt.freqrange), opt.freqrange = STUDY.etc.erspparams.freqrange; end;
end;
if ~isempty(opt.channels), fileExt = [ '.dat' tmpDataType ];
else                       fileExt = [ '.ica' tmpDataType ];
end;

% first subject data file
% -----------------------
testSubjectFile = fullfile(ALLEEG(1).filepath, [ ALLEEG(1).subject fileExt ]);

% list of subjects
% ----------------
allSubjects = { STUDY.datasetinfo.subject };
uniqueSubjects = unique(allSubjects);
STUDY.subject = uniqueSubjects;
if ischar(opt.subject) && ~isempty(opt.subject), subjectList = {opt.subject}; else subjectList = opt.subject; end;
if isempty(subjectList)
    if isnan(opt.design), subjectList = STUDY.subject;
    else subjectList = STUDY.design(opt.design).cases.value; 
    end;
end;

% options
% -------
if strcmpi(dtype, 'erp'), opts = { 'timelimits', opt.timerange };
else                      opts = { 'freqlimits', opt.freqrange };
end;
opts = { opts{:} 'singletrials' opt.singletrials };

for iSubj = 1:length(subjectList)
    
    % check cache
    bigstruct = [];
    if ~isempty(opt.channels), bigstruct.channel = opt.channels;
    else                       bigstruct.cluster = opt.clusters; % there can only be one cluster
    end;
    bigstruct.datatype     = opt.datatype;
    bigstruct.singletrials = opt.singletrials;
    bigstruct.subject      = subjectList{iSubj};
    bigstruct.component    = opt.component;
    bigstruct.subbaseline  = opt.subbaseline;
    bigstruct.options      = opts;
    if isnan(opt.design)
         bigstruct.design.variable = struct([]);
    else bigstruct.design.variable = STUDY.design(opt.design).variable;
    end;

    % find component indices
    % ----------------------
    if ~isempty(opt.clusters)
        datasetInds = strmatch(subjectList{iSubj}, { STUDY.datasetinfo.subject }, 'exact');
        compList    = [];
        for iDat = datasetInds(:)'
            indSet   = find(STUDY.cluster(opt.clusters).sets(1,:) == iDat); % each column contain info about the same subject
            compList = [ compList STUDY.cluster(opt.clusters).comps(indSet)' ]; % so we many only consider the first row
        end;
    end;
    
    % read all channels/components at once
    hashcode = gethashcode(std_serialize(bigstruct));
    [STUDY.cache tmpstruct] = eeg_cache(STUDY.cache, hashcode);

    if ~isempty(tmpstruct)
        dataTmp{iSubj}   = tmpstruct{1};
        xvals            = tmpstruct{2};
        yvals            = tmpstruct{3};
        eventsTmp{iSubj} = tmpstruct{4};
        params           = tmpstruct{5};
    else
        datInds = find(strncmp( subjectList{iSubj}, allSubjects, max(cellfun(@length, allSubjects))));
        fileName = fullfile(STUDY.datasetinfo(datInds(1)).filepath, [ subjectList{iSubj} fileExt ]);
        if ~isempty(opt.channels)
             [dataTmp{iSubj} params xvals yvals eventsTmp{iSubj} ] = std_readfile( fileName, 'designvar', bigstruct.design.variable, opts{:}, 'channels', opt.channels);
        else [dataTmp{iSubj} params xvals yvals eventsTmp{iSubj} ] = std_readfile( fileName, 'designvar', bigstruct.design.variable, opts{:}, 'components', compList);
        end;

        if ~strcmpi(opt.datatype, 'ersp') && ~strcmpi(opt.datatype, 'itc') % ERP or spectrum
            if strcmpi(opt.singletrials, 'off')
                dataTmp{iSubj} = cellfun(@(x)squeeze(mean(x,2)), dataTmp{iSubj}, 'uniformoutput', false);
            end;
        else
            dataTmp{iSubj} = cellfun(@(x)processtf(x, xvals, opt.datatype, opt.singletrials, params), dataTmp{iSubj}, 'uniformoutput', false);
        end;
        if ~isempty(eventsTmp{iSubj}{1}) && strcmpi(opt.singletrials, 'off')
            eventsTmp{iSubj} = cellfun(@(x)squeeze(mean(x)), eventsTmp{iSubj}, 'uniformoutput', false);
        end;
        STUDY.cache = eeg_cache(STUDY.cache, hashcode, { dataTmp{iSubj} xvals yvals eventsTmp{iSubj} params });
    end;
end;

% if single trials put channels in 2nd dim and trials in last dim
if strcmpi(opt.singletrials, 'on') && length(opt.channels) > 1
    for iCase = 1:length(dataTmp)
        for iItem = 1:length(dataTmp{1}(:))
            dataTmp{iCase}{iItem} = permute(dataTmp{iCase}{iItem}, [1 3 2]);
        end;
    end;
end;

% store data for all subjects
if strcmpi(opt.datatype, 'erp') || strcmpi(opt.datatype, 'spec')
     if length(opt.channels) > 1, dim = 3; else dim = 2; end;
else if length(opt.channels) > 1, dim = 4; else dim = 3; end;
end;
datavals = reorganizedata(dataTmp, dim);
events   = reorganizedata(eventsTmp, 1);

% fix component polarity if necessary
% -----------------------------------
componentPol = [];
if isempty(opt.channels) && strcmpi(dtype, 'erp') && isempty(opt.channels) && strcmpi(opt.componentpol, 'on')
    disp('Reading component scalp topo polarities - this is done to invert some ERP component polarities');
    STUDY = std_readtopoclust(STUDY, ALLEEG, opt.clusters);
    componentPol = STUDY.cluster(opt.clusters).topopol;
    if isempty(componentPol)
        disp('Cluster topographies absent - cannot adjust single component ERP polarities');
    end;
    for iItem = 1:length(datavals)
        datavals{iItem} = bsxfun(@times, datavals{iItem}, componentPol);
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

% reorganize data
% ---------------
function datavals = reorganizedata(dataTmp, dim)
    datavals = cell(size(dataTmp{1}));
    for iItem=1:length(dataTmp{1}(:)')
        numItems = sum(cellfun(@(x)size(x{iItem},dim), dataTmp));
        switch dim
            case 2, datavals{iItem} = zeros([ size(dataTmp{1}{iItem},1) numItems], 'single'); 
            case 3, datavals{iItem} = zeros([ size(dataTmp{1}{iItem},1) size(dataTmp{1}{iItem},2) numItems], 'single'); 
            case 4, datavals{iItem} = zeros([ size(dataTmp{1}{iItem},1) size(dataTmp{1}{iItem},2) size(dataTmp{1}{iItem},3) numItems], 'single'); 
        end;
    end;
    for iItem=1:length(dataTmp{1}(:)')
        count = 1;
        for iCase = 1:length(dataTmp)
            if ~isempty(dataTmp{iCase}{iItem})
                numItems = size(dataTmp{iCase}{iItem},dim);
                switch dim
                    case 2, datavals{iItem}(:,count:count+numItems-1) = dataTmp{iCase}{iItem}; 
                    case 3, datavals{iItem}(:,:,count:count+numItems-1) = dataTmp{iCase}{iItem}; 
                    case 4, datavals{iItem}(:,:,:,count:count+numItems-1) = dataTmp{iCase}{iItem};
                end;
                count = count+numItems;
            end;
        end;
    end;
    
% call newtimef (duplicate function in std_erspplot)
% --------------
function dataout = processtf(dataSubject, xvals, datatype, singletrials, g)

    % compute ITC or ERSP
    if strcmpi(datatype, 'ersp')
        P = dataSubject .* conj(dataSubject);
        dataout = newtimeftrialbaseln(P, xvals, g);
        if strcmpi(singletrials, 'off')
            if ndims(dataout) == 4,     dataout = mean(dataout, 4);
            elseif ndims(dataout) == 3, dataout = mean(dataout, 3);
            end;
        end;
    else
        if strcmpi(singletrials, 'off')
            if ~isfield(params, 'itctype'), params.itctype = 'phasecoher'; end;
            dataout = newtimefitc(dataSubject, g.itctype);
        end;
    end;
    