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
% nc = max(length(STUDY.design(opt.design).variable(1).value),1);
% ng = max(length(STUDY.design(opt.design).variable(2).value),1);

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
if ~isempty(opt.channels), fileExt = '.daterp';
else                       fileExt = '.icaerp';
end;

% first subject data file
% -----------------------
testSubjectFile = fullfile(ALLEEG(1).filepath, [ ALLEEG(1).subject fileExt ]);

for ind = 1:length(finalinds) % scan channels or components

    % list of subjects
    % ----------------
    allSubjects = { STUDY.datasetinfo.subject };
    uniqueSubjects = unique(allSubjects);
    STUDY.subject = uniqueSubjects;
        
    % find indices
    % ------------
    if ~isempty(opt.channels)
        tmpstruct = STUDY.changrp(finalinds(ind));
    else
        tmpstruct = STUDY.cluster(finalinds(ind));
    end;

    % check if data is already here using hashcode
    % --------------------------------------------
    dataread = false;
    bigstruct = [];
    bigstruct.timerange = opt.timerange;
    bigstruct.freqrange = opt.freqrange;
    bigstruct.rmsubjmean   = opt.rmsubjmean;
    bigstruct.singletrials = opt.singletrials;
    bigstruct.subject      = opt.subject;
    bigstruct.design       = STUDY.design(opt.design);
    hashcode = gethashcode(std_serialize(bigstruct));
    if isfield(tmpstruct, [ dtype 'hashcode' ]) && strcmpi( tmpstruct.([ dtype 'hashcode' ]), hashcode), dataread = true; end;
    
    if dataread == false
        tmpstruct.([ dtype 'hashcode' ]) = hashcode;
        
        % reading options
        % ---------------
        fprintf([ 'Reading ' dtype ' data...\n' ]);
        if strcmpi(dtype, 'erp'), opts = { 'timelimits', opt.timerange };
        else                      opts = { 'freqlimits', opt.freqrange };
        end;
        
        % read the data and select channels
        % ---------------------------------
        subjectList = opt.subject;
        if isempty(subjectList), subjectList = STUDY.subject; end;
        for iSubj = length(subjectList):-1:1
            inds = find(strncmp( subjectList{iSubj}, allSubjects, max(cellfun(@length, allSubjects))));
            fileName = fullfile(STUDY.datasetinfo(inds(1)).filepath, [ subjectList{iSubj} fileExt ]); 
            [dataSubject{ iSubj } params xvals tmp events{ iSubj } ] = std_readfile( fileName,  'designvar', STUDY.design(opt.design).variable, opts{:}, 'channels', opt.channels(ind));
            % DEAL WITH COMPONENTS HERE
        end;
        
        % concatenate subject data
        if strcmpi(opt.singletrials, 'off')
            for iSubj = length(subjectList):-1:1
                for iCell = 1:length(dataSubject{1}(:))
                    alldata{  iCell}(:,iSubj) = mean(dataSubject{ iSubj }{ iCell },2);
                    if ~isempty(events{iSubj}{iCell})
                         allevents{iCell}(:,iSubj) = mean(events{      iSubj }{ iCell },2);
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
                         allevents{iCell}(:,allcounts(iCell)+1:allcounts(iCell)+cols) = events{      iSubj }{ iCell };
                    else allevents{iCell} = [];
                    end;
                    allcounts(iCell) = allcounts(iCell)+cols;
                end;
            end;
        end;
        alldata   = reshape(alldata  , size(dataSubject{1}));
        allevents = reshape(allevents, size(events{1}));
        
        % inverting component polaritites - HAVE TO CHECK HERE ABOUT THE NEW FRAMEWORK
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
        
        tmpstruct = setfield( tmpstruct, [ dtype 'data' ], alldata);
        tmpstruct = setfield( tmpstruct, [ dtype 'vars' ], allevents);
        if strcmpi(dtype, 'spec'), tmpstruct.specfreqs = xvals;
        else                       tmpstruct.erptimes  = xvals;
        end;
        
        % copy results to structure
        % -------------------------
        fieldnames = { [ dtype 'data' ] [ dtype 'vars' ] [ dtype 'freqs' ] [ dtype 'times' ] [ dtype 'hashcode' ] };
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
    datavals  = [];
    if strcmpi( dtype, 'spec'), xvals = getfield(structdat(allinds(1)), [ dtype 'freqs' ]);
    else                        xvals = getfield(structdat(allinds(1)), [ dtype 'times' ]);
    end;
   
    for chan = length(finalinds):-1:1
        tmpdat = getfield(structdat(finalinds(chan)), [ dtype 'data' ]); % only works for interpolated data
        for ind =  1:length(tmpdat(:))
            datavals{ind}(:,:,chan) = tmpdat{ind};
        end;
    end;
    datavals = reshape(datavals, size(tmpdat));
    
    for ind =  1:length(datavals(:))
        datavals{ind} = squeeze(permute(datavals{ind}, [1 3 2])); % time elec subjects
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


   