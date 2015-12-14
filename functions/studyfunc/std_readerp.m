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

function [STUDY, datavals, xvals] = std_readerp(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readerp;
    return;
end

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

dtype = opt.datatype;

% find channel indices
% --------------------
if ~isempty(opt.channels)
     allChangrp = lower({ STUDY.changrp.name });
     finalinds = std_chaninds(STUDY, opt.channels);
else finalinds = opt.clusters;
end;

% get the file extension
% ----------------------
if ~isempty(opt.channels), fileExt = [ '.dat' opt.datatype ];
else                       fileExt = [ '.ica' opt.datatype ];
end;

% first subject data file
% -----------------------
testSubjectFile = fullfile(ALLEEG(1).filepath, [ ALLEEG(1).subject fileExt ]);

newstruct = [];
for ind = 1:length(finalinds) % scan channels or clusters

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
    bigstruct.rmsubjmean   = opt.rmsubjmean;
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
        
        % reading options
        % ---------------
        fprintf([ 'Reading ' dtype ' data...\n' ]);
        if strcmpi(dtype, 'erp'), opts = { 'timelimits', opt.timerange };
        else                      opts = { 'freqlimits', opt.freqrange };
        end;
        
        % get component polarity if necessary
        % -----------------------------------
        componentPol = [];
        if isempty(opt.channels) && strcmpi(dtype, 'erp') && isempty(opt.channels) && strcmpi(opt.componentpol, 'on')
            disp('Reading component scalp topo polarities - this is done to invert some ERP component polarities');
            STUDY = std_readtopoclust(STUDY, ALLEEG, finalinds(ind));
            componentPol = STUDY.cluster(finalinds(ind)).topopol;
            if isempty(componentPol)
                disp('Cluster topographies absent - cannot adjust single component ERP polarities');
            end;
        end;
        
        % read the data and select channels
        % ---------------------------------
        if ischar(opt.subject) && ~isempty(opt.subject), subjectList = {opt.subject}; else subjectList = opt.subject; end;
        if isempty(subjectList), subjectList = bigstruct.design.cases.value; end;
        count = 1;
        for iSubj = length(subjectList):-1:1
            datInds = find(strncmp( subjectList{iSubj}, allSubjects, max(cellfun(@length, allSubjects))));
            fileName = fullfile(STUDY.datasetinfo(datInds(1)).filepath, [ subjectList{iSubj} fileExt ]); 
            
            if ~isempty(opt.channels)
                [dataSubject{ iSubj } params xvals tmp events{ iSubj } ] = std_readfile( fileName, 'designvar', bigstruct.design.variable, opts{:}, 'channels', opt.channels(ind));
            else
                % find components for a given cluster and subject
                setInds = [];
                for iDat = 1:length(datInds), setInds = [setInds find(STUDY.cluster(finalinds(ind)).sets(1,:) == datInds(iDat))' ]; end;
                if ~isempty(opt.component), setInds = intersect(setInds, opt.component); end;
                for iComp = 1:length(setInds)
                    comps = STUDY.cluster(finalinds(ind)).comps( setInds(iComp) );
                    [dataSubject{ count } params xvals tmp events{ count } ] = std_readfile( fileName,  'designvar', bigstruct.design.variable, opts{:}, 'components', comps);
                    if ~isempty(componentPol), for iCell = 1:length(dataSubject{ count }(:)), dataSubject{ count }{ iCell } = dataSubject{ count }{ iCell }*componentPol(setInds(iComp)); end; end;
                    count = count+1;
                end;
            end;
        end;

        % concatenate data - compute average if not dealing with (processing) single trials
        % ---------------------------------------------------------------------------------
        if strcmpi(opt.singletrials, 'off')
            for iSubj = length(dataSubject(:)):-1:1
                for iCell = 1:length(dataSubject{1}(:))
                    if isempty(dataSubject{ iSubj }{ iCell })
                        error(sprintf('Subject %s missing one experimental condition, remove subject and try again'));
                    end;
                    alldata{  iCell}(:,iSubj) = mean(dataSubject{ iSubj }{ iCell },2);
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
                    alldim(iCell) = alldim(iCell)+size(dataSubject{ iSubj }{ iCell },2);
                end;
            end;
            % initialize arrays
            for iCell = 1:length(dataSubject{1}(:))
                alldata{  iCell} = zeros(size(dataSubject{ 1 }{ 1 },1), alldim(iCell));
                allevents{iCell} = zeros(size(events{      1 }{ 1 },1), alldim(iCell));
            end;
            % populate arrays
            allcounts = zeros(size(dataSubject{1}));
            for iSubj = length(subjectList):-1:1
                for iCell = 1:length(dataSubject{1}(:))
                    cols = size(dataSubject{ iSubj }{ iCell },2);
                    alldata{  iCell}(:,allcounts(iCell)+1:allcounts(iCell)+cols) = dataSubject{ iSubj }{ iCell };
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

        % remove mean of each subject across groups and conditions - HAVE TO CHECK HERE ABOUT THE NEW FRAMEWORK
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
        
        newstruct(ind).([ dtype 'data' ]) = alldata;
        newstruct(ind).([ dtype 'vars' ]) = allevents;
        if strcmpi(dtype, 'spec'), newstruct(ind).specfreqs = xvals;
        else                       newstruct(ind).erptimes  = xvals;
        end;
        STUDY.cache = eeg_cache(STUDY.cache, hashcode, newstruct(ind));
        
    end;
end;

% if several channels, agregate them
% ----------------------------------
allinds   = finalinds;
if strcmpi( dtype, 'spec'), xvals = newstruct(1).([ dtype 'freqs' ]);
else                        xvals = newstruct(1).([ dtype 'times' ]);
end;
if ~isempty(opt.channels)
    % concatenate channels if necessary
    datavals  = [];   
    for chan = length(finalinds):-1:1
        tmpdat = newstruct(chan).([ dtype 'data' ]); % only works for interpolated data
        for ind =  1:length(tmpdat(:))
            datavals{ind}(:,:,chan) = tmpdat{ind};
        end;
    end;
    datavals = reshape(datavals, size(tmpdat));
    
    for ind =  1:length(datavals(:))
        datavals{ind} = squeeze(permute(datavals{ind}, [1 3 2])); % time elec subjects
    end;
else
    % in practice clusters are read one by one
    % so it is fine to take the first element only
    datavals = newstruct(1).([ dtype 'data' ]);
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


   