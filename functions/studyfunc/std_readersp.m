% std_readersp() - load ERSP measures for data channels or  for all 
%                  components of a specified cluster. This function is 
%                  also being used to read ITC and ERPimage data.
% Usage:
%   >> [STUDY, erspdata, times, freqs, erspbase] = ...
%                   std_readersp(STUDY, ALLEEG, varargin);
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
%  'forceread' - ['on'|'off'] Force rereading data from disk.
%                Default is 'off'.
%  'subject'   - [string] select a specific subject {default:all}
%  'component' - [integer] select a specific component in a cluster.
%                This is the index of the component in the cluster not the
%                component number {default:all}
%  'datatype'  - {'ersp'|'itc'|'erpim'} This function is used to read all 
%                2-D STUDY matrices stored on disk (not only ERSP). It may
%                read ERSP ('ersp' option), ITC ('itc' option) or ERPimage
%                data ('erpim' option).
%
% ERSP specific options:
%  'timerange' - [min max] time range {default: whole measure range}
%  'freqrange' - [min max] frequency range {default: whole measure range}
%  'subbaseline' - ['on'|'off'] subtract the ERSP baseline for paired
%                conditions. The conditions for which baseline is removed
%                are indicated on the command line. See help
%                std_studydesign for more information about paired and
%                unpaired variables.
%
% ERPimage specific option:
%   This function is used to read all 2-D STUDY matrices stored on disk
%   (this includes ERPimages). It therefore takes as input specific 
%   ERPimage options. Note that the 'singletrials' optional input is
%   irrelevant for ERPimages (which are always stored as single trials).
%   'concatenate' - ['on'|'off'] read concatenated ERPimage data ('on') or 
%                   stacked ERPimage data. See help std_erpimage for more 
%                   information.
%   'timerange'   - [min max] time range {default: whole measure range}
%   'trialrange'  - [min max] read only a specific range of the ERPimage
%                   output trials {default: whole measure range}
% Output:
%  STUDY    - updated studyset structure
%  erspdata - [cell array] ERSP data (the cell array size is 
%             condition x groups). This may also be ITC data or ERPimage
%             data (see above).
%  times    - [float array] array of time points
%  freqs    - [float array] array of frequencies. For ERPimage this
%             contains trial indices.
%  erspbase - [cell array] baseline values
%  events   - [cell array] events (ERPimage only).
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

% TO INVESTIGAGE
% TO INVESTIGAGE
% TO INVESTIGAGE
% TO INVESTIGAGE
% TO INVESTIGAGE
% TO INVESTIGAGE
% TO INVESTIGAGE

% REMOVE COMMON BASELINE AS A REALTIME OPTION

% TO INVESTIGAGE
% TO INVESTIGAGE
% TO INVESTIGAGE
% TO INVESTIGAGE
% TO INVESTIGAGE
% TO INVESTIGAGE

function [STUDY, erspdata, alltimes, allfreqs, erspbase, events, unitPower] = std_readersp(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readersp;
    return;
end
events = {};
unitPower = 'dB';
erspbase = [];

STUDY = pop_erspparams( STUDY, 'default');
STUDY = pop_erpimparams(STUDY, 'default');
[opt moreopts] = finputcheck( varargin, { ...
    'design'        'integer' []             STUDY.currentdesign;    
    'channels'      'cell'    []             {};
    'clusters'      'integer' []             [];
    'trialrange'    'real'    []             STUDY.etc.erpimparams.trialrange;
    'freqrange'     'real'    []             STUDY.etc.erspparams.freqrange;
    'timerange'     'real'    []             NaN;
    'singletrials'  'string'  { 'on','off' } 'off';
    'forceread'     'string'  { 'on','off' } 'off';
    'concatenate'   'string'  { 'on','off' } STUDY.etc.erpimparams.concatenate;
    'subbaseline'   'string'  { 'on','off' } STUDY.etc.erspparams.subbaseline;
    'component'     'integer' []             [];
    'infotype'      'string'  { 'ersp','itc','erpim' } 'ersp'; ... % deprecated
    'datatype'      'string'  { 'ersp','itc','erpim' } 'ersp'; ...
    'subject'       'string'  []             '' }, ...
    'std_readersp', 'ignore');
if isstr(opt), error(opt); end;
if ~strcmpi(opt.infotype, 'ersp'), opt.datatype = opt.infotype; end;
if strcmpi(opt.datatype, 'erpim'), if isnan(opt.timerange), opt.timerange = STUDY.etc.erpimparams.timerange; end;
                                   opt.freqrange = opt.trialrange;
                                   ordinate      = 'trials';
else                               if isnan(opt.timerange) opt.timerange = STUDY.etc.erspparams.timerange; end;
                                   ordinate      = 'freqs';
end;
if ~isnan(opt.design) && length(STUDY.design(opt.design).variable) > 0, paired1 = STUDY.design(opt.design).variable(1).pairing; else paired1 = 'off'; end;
if ~isnan(opt.design) && length(STUDY.design(opt.design).variable) > 1, paired2 = STUDY.design(opt.design).variable(2).pairing; else paired2 = 'off'; end;
dtype = opt.datatype;

% find channel indices
% --------------------
if ~isempty(opt.channels)
     allChangrp = lower({ STUDY.changrp.name });
     finalinds = std_chaninds(STUDY, opt.channels);
else finalinds = opt.clusters;
end;

newstruct = [];
for ind = 1:length(finalinds)
    %erspbase = cell( nc, ng );

    % list of subjects
    % ----------------
    allSubjects = { STUDY.datasetinfo.subject };
    uniqueSubjects = unique(allSubjects);
    STUDY.subject = uniqueSubjects;

    % check if data is already here using hashcode
    % --------------------------------------------
    bigstruct = [];
    if ~isempty(opt.channels), bigstruct.channel = opt.channels{ind};
    else                       bigstruct.cluster = opt.clusters(ind);
    end;
    bigstruct.datatype     = opt.datatype;
    bigstruct.timerange    = opt.timerange;
    bigstruct.freqrange    = opt.freqrange;
    bigstruct.trialrange   = opt.trialrange;
    %bigstruct.rmsubjmean   = opt.rmsubjmean;
    bigstruct.singletrials = opt.singletrials;
    bigstruct.subject      = opt.subject;
    bigstruct.component    = opt.component;
    if isnan(opt.design)
         bigstruct.design.cases.value = STUDY.subject;
         bigstruct.design.variable    = struct([]);
    else bigstruct.design  = STUDY.design(opt.design);
    end;
    hashcode = gethashcode(std_serialize(bigstruct));
    [STUDY.cache tmpstruct] = eeg_cache(STUDY.cache, hashcode);
    
    if ~isempty(tmpstruct)
        if isempty(newstruct), newstruct = tmpstruct;
        else                   newstruct(ind) = tmpstruct;
        end;
    else
        % reserve arrays
        % --------------
        events   = {};
        %ersp     = cell( nc, ng );
        %erspinds = cell( nc, ng );
        
        % reading options
        % ---------------
        fprintf([ 'Reading ' dtype ' data...\n' ]);
        if ~strcmpi(dtype, 'erpim')
            opts = { 'measure', 'timef' 'freqlimits', opt.freqrange }; % 'timelimits', opt.timerange, (time is selected later to allow for baseline removal)
        else
            opts = { 'measure', 'erpim' 'triallimits', opt.trialrange 'timerange' opt.timerange };
        end;
        
        % read the data and select channels
        % ---------------------------------
        if ischar(opt.subject) && ~isempty(opt.subject), subjectList = {opt.subject}; else subjectList = opt.subject; end;
        if isempty(subjectList), subjectList = STUDY.design(STUDY.currentdesign).cases.value; end;
        count = 1;
        for iSubj = 1:length(subjectList)
            datInds = find(strncmpi( subjectList{iSubj}, allSubjects, max(cellfun(@length, allSubjects))));
            
            if ~isempty(opt.channels)
                fileName = fullfile(STUDY.datasetinfo(datInds(1)).filepath, [ subjectList{iSubj} '.dat' fastif(strcmpi(dtype, 'erpim'), 'erpim', 'timef') ]);
                [dataSubject{ iSubj } params xvals yvals events{ iSubj } ] = std_readfile( fileName,  'designvar', bigstruct.design.variable, opts{:}, 'channels', opt.channels(ind));
            else
                % find components for a given cluster and subject
                fileName = fullfile(STUDY.datasetinfo(datInds(1)).filepath, [ subjectList{iSubj} '.ica' fastif(strcmpi(dtype, 'erpim'), 'erpim', 'timef') ]);
                setInds = [];
                for iDat = 1:length(datInds), setInds = [setInds find(STUDY.cluster(finalinds(ind)).sets(1,:) == datInds(iDat))' ]; end;
                if ~isempty(opt.component), setInds = intersect(setInds, opt.component); end;
                for iComp = 1:length(setInds)
                    comps = STUDY.cluster(finalinds(ind)).comps( setInds(iComp) );
                    [dataSubject{ count } params xvals yvals events{ count } ] = std_readfile( fileName,  'designvar', bigstruct.design.variable, opts{:}, 'components', comps);
                    count = count+1;
                end;
            end;
        end;
        
        % Issues with ERSP
        
        % when reading the data, if we select a specific frequency range,
        % this is going to change the baseline which is computed from
        % the available visible baseline.
        
        
        % concatenate data - compute average if not dealing with (processing) single trials
        % ---------------------------------------------------------------------------------
        if strcmpi(opt.singletrials, 'off')
            for iSubj = length(dataSubject(:)):-1:1
                for iCell = 1:length(dataSubject{1}(:))
                    if isempty(dataSubject{ iSubj }{ iCell })
                        error(sprintf('Subject %s missing one experimental condition, remove subject and try again'));
                    end;
                    tmpdat = callnewtimef(dataSubject{ iSubj }{ iCell }, xvals, yvals, ALLEEG(1).pnts, [ALLEEG(1).xmin ALLEEG(1).xmax]*1000, opt.datatype, params);
                    alldata{  iCell}(:,:,iSubj) = tmpdat;
                    
                    if ~isempty(events{iSubj}{iCell})
                         allevents{iCell}(:,iSubj) = mean(events{ iSubj }{ iCell },2);
                    else allevents{iCell} = [];
                    end;
                end;
            end;
        else
            % calculate dimensions
            alldim = zeros(size(dataSubject{1}));
            for iSubj = length(subjectList):-1:1
                for iCell = 1:length(dataSubject{1}(:))
                    alldim(iCell) = alldim(iCell)+size(dataSubject{ iSubj }{ iCell },3);
                end;
            end;
            % initialize arrays
            for iCell = 1:length(dataSubject{1}(:))
                alldata{  iCell} = zeros(size(dataSubject{ 1 }{ 1 },1), size(dataSubject{ 1 }{ 1 },2), alldim(iCell));
                allevents{iCell} = zeros(size(events{      1 }{ 1 },1), alldim(iCell));
            end;
            % populate arrays
            allcounts = zeros(size(dataSubject{1}));
            for iSubj = length(subjectList):-1:1
                for iCell = 1:length(dataSubject{1}(:))
                    cols = size(dataSubject{ iSubj }{ iCell },3);
                    alldata{  iCell}(:,:,allcounts(iCell)+1:allcounts(iCell)+cols) = dataSubject{ iSubj }{ iCell };
                    if ~isempty(events{iSubj}{iCell})
                        allevents{iCell}(:,allcounts(iCell)+1:allcounts(iCell)+cols) = events{ iSubj }{ iCell };
                    else allevents{iCell} = [];
                    end;
                    allcounts(iCell) = allcounts(iCell)+cols;
                end;
            end;
        end;
        alldata   = reshape(alldata  , size(dataSubject{1}));
        allevents = reshape(allevents, size(events{1}));
        
        newstruct(ind).([ dtype 'data' ])   = alldata;
        newstruct(ind).([ dtype 'vars' ])   = allevents;
        newstruct(ind).([ dtype 'times' ])  = xvals;
        newstruct(ind).([ dtype 'freqs' ])  = yvals;
        newstruct(ind).([ dtype 'params' ]) = params;
        STUDY.cache = eeg_cache(STUDY.cache, hashcode, newstruct(ind));
        
    end;
end;

% output units
% -----------
tmpparams = newstruct(1).([ dtype 'params' ]);
if ~isfield(tmpparams, 'baseline'), tmpparams.baseline = 0;     end;
if ~isfield(tmpparams, 'scale'   ), tmpparams.scale    = 'log'; end;
if ~isfield(tmpparams, 'basenorm'), tmpparams.basenorm = 'off'; end;
if strcmpi(tmpparams.scale, 'log')
    if strcmpi(tmpparams.basenorm, 'on')
        unitPower = '10*log(std.)'; % impossible
    elseif isnan(tmpparams.baseline)
        unitPower = '10*log10(\muV^{2}/Hz)';
    else
        unitPower = 'dB';
    end;
else
    if strcmpi(tmpparams.basenorm, 'on')
        unitPower = 'std.';
    elseif isnan(tmpparams.baseline)
        unitPower = '\muV^{2}/Hz';
    else
        unitPower = '% of baseline';
    end;
end;

% return structure
% ----------------
if nargout <2
    return
end
allinds   = finalinds;
allfreqs = newstruct(1).([ dtype 'freqs']);
alltimes = newstruct(1).([ dtype 'times']);
if ~isempty(opt.channels)
        % concatenate channels if necessary
    erspdata  = [];
    for chan = length(finalinds):-1:1
        tmpdat = newstruct(chan).([ dtype 'data' ]); % only works for interpolated data
        for ind =  1:length(tmpdat(:))
            erspdata{ind}(:,:,:,chan) = tmpdat{ind};
        end;
    end;
    erspdata = reshape(erspdata, size(tmpdat));
    
    for ind =  1:length(erspdata(:))
        erspdata{ind} = squeeze(permute(erspdata{ind}, [1 2 4 3])); % time freq elec subjects
    end;
else
    % in practice clusters are read one by one
    % so it is fine to take the first element only
    erspdata = newstruct(1).([ dtype 'data' ]);
end;

% get events
events = newstruct(1).([ dtype 'vars' ]); % only works for interpolated data

% compute ERSP baseline
% ---------------------
function meanpowbase = computeerspbaseline(erspbase, singletrials)

    try
        len = length(erspbase(:));
        count = 0;
        for index = 1:len
            if ~isempty(erspbase{index})
                if strcmpi(singletrials, 'on')
                    if count == 0, meanpowbase = abs(mean(erspbase{index},3));
                    else           meanpowbase = meanpowbase + abs(mean(erspbase{index},3));
                    end;
                else
                    if count == 0, meanpowbase = abs(erspbase{index});
                    else           meanpowbase = meanpowbase + abs(erspbase{index});
                    end;
                end;
                count = count+1;
            end;
        end;
        meanpowbase = reshape(meanpowbase  , [size(meanpowbase,1) 1 size(meanpowbase,2)])/count;
    catch,
        error([ 'Problem while subtracting common ERSP baseline.' 10 ...
                'Common baseline subtraction is performed based on' 10 ...
                'pairing settings in your design. Most likelly, one' 10 ...
                'independent variable should not have its data paired.' ]);
    end;
        
% remove ERSP baseline
% ---------------------
function ersp = removeerspbaseline(ersp, erspbase, meanpowbase)
    convert2log = 0;
    for g = 1:size(ersp,2)        % ng = number of groups
        for c = 1:size(ersp,1)
            if ~isempty(erspbase{c,g}) && ~isempty(ersp{c,g})
                erspbasetmp = reshape(erspbase{c,g}, [size(meanpowbase,1) 1 size(meanpowbase,3)]);
                if any(erspbasetmp(:) > 1000)
                    convert2log = 1;
                end;
                tmpmeanpowbase = repmat(meanpowbase, [1 size(ersp{c,g},2) 1]);
                if convert2log
                     ersp{c,g} = ersp{c,g} - repmat(10*log10(erspbasetmp), [1 size(ersp{c,g},2) 1 1]) + 10*log10(tmpmeanpowbase);
                else ersp{c,g} = ersp{c,g} - repmat(abs(erspbasetmp), [1 size(ersp{c,g},2) 1 1]) + tmpmeanpowbase;
                end;
            end;
        end;
    end;

% call newtimef (duplicate function in std_erspplot)
% --------------
function [dataout tmpParams] = callnewtimef(dataSubject, xvals, yvals, pnts, tlimits, datatype, params);

    precomputed.tfdata = dataSubject;
    precomputed.times = xvals;
    precomputed.freqs = yvals;
    precomputed.recompute = datatype;

    cycles = params.cycles;
    params = rmfield(params, 'cycles');
    tmpParams = fieldnames(params)';
    tmpParams(2,:) = struct2cell(params)';
    srate = 1;
    
    if strcmpi(datatype, 'ersp')
        dataout = newtimef([],pnts,tlimits, srate, cycles, 'precomputed', precomputed, 'verbose', 'off', tmpParams{:});
    else
        [tmp dataout] = newtimef([],pnts,tlimits, srate, cycles, 'precomputed', precomputed, 'verbose', 'off', tmpParams{:});
    end;
