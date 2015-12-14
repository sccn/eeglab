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
        % reading options
        % ---------------
        fprintf([ 'Reading ' dtype ' data...\n' ]);
        if strcmpi(dtype, 'erp'), opts = { 'timelimits', opt.timerange };
        else                      opts = { 'freqlimits', opt.freqrange };
        end;
        
        if strcmpi(dtype, 'erpim') % ERP images
            
            if strcmpi(opt.singletrials, 'on')
                error( [ 'Single trial loading option not supported with STUDY ERP-image' 10 '(there is no such thing as a single-trial ERPimage)' ]);
            end;
            % read the data and select channels
            % ---------------------------------
            setinfo = STUDY.design(opt.design).cell;
            erpim    = cell( nc, ng );
            events   = cell( nc, ng );
            for c = 1:nc
                for g = 1:ng
                    if ~isempty(setinds{c,g})
                        if ~isempty(opt.channels), opts = { 'channels',    allChangrp(allinds{c,g}(:)), 'timelimits', opt.timerange, 'triallimits', opt.trialrange, 'concatenate', opt.concatenate };
                        else                       opts = { 'components',  allinds{c,g}(:), 'timelimits', opt.timerange, 'triallimits', opt.trialrange, 'concatenate', opt.concatenate };
                        end;
                        [erpim{c, g} tmpparams alltimes alltrials events{c,g}] = std_readfile( setinfo(setinds{c,g}(:)), 'measure', 'erpim', opts{:});
                        fprintf('.');
                    end;
                end;
            end;
            if strcmpi(opt.concatenate, 'on'), alltrials = []; end;
            fprintf('\n');
            ersp     = erpim;
            allfreqs = alltrials;
            
        else % ERSP/ITC
            
            % reserve arrays
            % --------------
            events   = {};
            %ersp     = cell( nc, ng );
            %erspinds = cell( nc, ng );
            
            % reading options
            % ---------------
            fprintf([ 'Reading ' dtype ' data...\n' ]);
            opts = { 'measure', 'timef' 'freqlimits', opt.freqrange }; % 'timelimits', opt.timerange, (time is selected later to allow for baseline removal)
            
            % read the data and select channels
            % ---------------------------------
            if ischar(opt.subject) && ~isempty(opt.subject), subjectList = {opt.subject}; else subjectList = opt.subject; end;
            if isempty(subjectList), subjectList = STUDY.design(STUDY.currentdesign).cases.value; end;
            count = 1;
            for iSubj = length(subjectList):-1:1
                datInds = find(strncmp( subjectList{iSubj}, allSubjects, max(cellfun(@length, allSubjects))));
                fileName = fullfile(STUDY.datasetinfo(datInds(1)).filepath, [ subjectList{iSubj} '.dattimef' ]);
                
                if ~isempty(opt.channels)
                    [dataSubject{ iSubj } params xvals yvals events{ iSubj } ] = std_readfile( fileName,  'designvar', bigstruct.design.variable, opts{:}, 'channels', opt.channels(ind));
                else
                    % find components for a given cluster and subject
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
            
            % Issues
            
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
                        tmpdat = callnewtimef(dataSubject{ iSubj }{ iCell }, xvals, yvals, ALLEEG(1).pnts, [ALLEEG(1).xmin ALLEEG(1).xmax]*1000, opt.datatype, opt.timerange, params);
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
            
            newstruct(ind).([ dtype 'data' ])   = alldata;
            newstruct(ind).([ dtype 'vars' ])   = allevents;
            newstruct(ind).([ dtype 'times' ])  = xvals;
            newstruct(ind).([ dtype 'freqs' ])  = yvals;
            newstruct(ind).([ dtype 'params' ]) = params;
            STUDY.cache = eeg_cache(STUDY.cache, hashcode, newstruct(ind));
        
%            alldata   = reshape(alldata  , size(dataSubject{1}));
%            allevents = reshape(allevents, size(events{1}));
        
%             % read the data and select channels
%             % ---------------------------------
%             fprintf('Reading all %s data:', upper(dtype));
%             setinfo = STUDY.design(opt.design).cell;
%             if strcmpi(opt.singletrials, 'on')
%                 for c = 1:nc
%                     for g = 1:ng
%                         if ~isempty(opt.channels)
%                             allsubjects = { STUDY.design(opt.design).cell.case };
%                             if ~isempty(opt.subject), inds = strmatch( opt.subject, allsubjects(setinds{c,g}));
%                             else inds = 1:length(allinds{c,g}); end;
%                         else
%                             if ~isempty(opt.component) inds = find( allinds{c,g} == STUDY.cluster(finalinds(ind)).comps(opt.component));
%                             else inds = 1:length(allinds{c,g}); end;
%                         end;
%                         if ~isempty(inds)
%                             count{c, g} = 1;
%                             for indtmp = 1:length(inds)
%                                 setindtmp = STUDY.design(opt.design).cell(setinds{c,g}(inds(indtmp)));
%                                 tmpopts = { 'measure', 'timef' 'timelimits', opt.timerange, 'freqlimits', opt.freqrange };
%                                 if ~isempty(opt.channels), [ tmpersp tmpparams alltimes allfreqs] = std_readfile(setindtmp, 'channels', allChangrp(allinds{c,g}(inds(indtmp))), tmpopts{:});
%                                 else                       [ tmpersp tmpparams alltimes allfreqs] = std_readfile(setindtmp, 'components',          allinds{c,g}(inds(indtmp)),  tmpopts{:});
%                                 end;
%                                 indices = [count{c, g}:count{c, g}+size(tmpersp,3)-1];
%                                 if indtmp == 1
%                                     ersp{c, g} = permute(tmpersp, [2 1 3]);
%                                 else ersp{c, g}(:,:,indices) = permute(tmpersp, [2 1 3]);
%                                 end;
%                                 erspinds{c, g}(1:2,indtmp) = [ count{c, g} count{c, g}+size(tmpersp,3)-1 ];
%                                 count{c, g} = count{c, g}+size(tmpersp,3);
%                                 if size(tmpersp,3) ~= sum(cellfun(@length, setindtmp.trials))
%                                     error( sprintf('Wrong number of trials in datafile for design index %d\n', setinds{c,g}(inds(indtmp))));
%                                 end;
%                             end;
%                         end;
%                     end;
%                 end;
%             else
%                 for c = 1:nc
%                     for g = 1:ng
%                         if ~isempty(setinds{c,g})
%                             if ~isempty(opt.channels), opts = { 'channels',  allChangrp(allinds{c,g}(:)), 'timelimits', opt.timerange, 'freqlimits', opt.freqrange };
%                             else                       opts = { 'components',           allinds{c,g}(:) , 'timelimits', opt.timerange, 'freqlimits', opt.freqrange };
%                             end;
%                             if strcmpi(dtype, 'ersp')
%                                  erspbase{c, g}                             = std_readfile( setinfo(setinds{c,g}(:)), 'measure', 'erspbase', opts{:});
%                                  [ ersp{c, g} tmpparams alltimes allfreqs ] = std_readfile( setinfo(setinds{c,g}(:)), 'measure', 'ersp'    , opts{:});
%                             else [ ersp{c, g} tmpparams alltimes allfreqs ] = std_readfile( setinfo(setinds{c,g}(:)), 'measure', 'itc'     , opts{:});
%                                  ersp{c, g} = abs(ersp{c, g});
%                             end;
%                             fprintf('.');
%                             %ersp{c, g}      = permute(ersp{c, g}           , [3 2 1]);
%                             %erspbase{c, g}  = 10*log(permute(erspbase{c, g}, [3 2 1]));
%                         end;
%                     end;
%                 end;
%             end;
%             fprintf('\n');
%             
%             % compute ERSP or ITC if trial mode
%             % (since only the timef have been loaded)
%             % ---------------------------------------
%             if strcmpi(opt.singletrials, 'on')
%                 tmpparams2      = fieldnames(tmpparams); tmpparams2 = tmpparams2';
%                 tmpparams2(2,:) = struct2cell(tmpparams);
%                 precomp.times = alltimes;
%                 precomp.freqs = allfreqs;
%                 precomp.recompute = dtype;
%                 for c = 1:nc
%                     for g = 1:ng
%                         if ~isempty(ersp{c,g})
%                             precomp.tfdata = permute(ersp{c,g}, [2 1 3]);
%                             if strcmpi(dtype, 'itc')
%                                 [tmp ersp{c,g}] = newtimef(zeros(ALLEEG(1).pnts,2), ALLEEG(1).pnts, [ALLEEG(1).xmin ALLEEG(1).xmax]*1000, ...
%                                     ALLEEG(1).srate, [], tmpparams2{:}, 'precomputed', precomp, 'verbose', 'off');
%                             elseif strcmpi(dtype, 'ersp')
%                                 [ersp{c,g} tmp] = newtimef(zeros(ALLEEG(1).pnts,2), ALLEEG(1).pnts, [ALLEEG(1).xmin ALLEEG(1).xmax]*1000, ...
%                                     ALLEEG(1).srate, [], tmpparams2{:}, 'precomputed', precomp, 'verbose', 'off');
%                             end;
%                             ersp{c,g} = permute(ersp{c,g}, [2 1 3]);
%                         end;
%                     end;
%                 end;
%             end;
%             
%             % compute average baseline across groups and conditions
%             % -----------------------------------------------------
%             if strcmpi(opt.subbaseline, 'on') && strcmpi(dtype, 'ersp')
%                 if strcmpi(opt.singletrials, 'on')
%                     disp('WARNING: no ERSP baseline may not be subtracted when using single trials');
%                 else
%                     disp('Recomputing baseline...');
%                     if strcmpi(paired1, 'on') && strcmpi(paired2, 'on')
%                         disp('Removing ERSP baseline for both indep. variables');
%                         meanpowbase = computeerspbaseline(erspbase(:), opt.singletrials);
%                         ersp        = removeerspbaseline(ersp, erspbase, meanpowbase);
%                     elseif strcmpi(paired1, 'on')
%                         disp('Removing ERSP baseline for first indep. variables (second indep. var. is unpaired)');
%                         for g = 1:ng        % ng = number of groups
%                             meanpowbase = computeerspbaseline(erspbase(:,g), opt.singletrials);
%                             ersp(:,g)   = removeerspbaseline(ersp(:,g), erspbase(:,g), meanpowbase);
%                         end;
%                     elseif strcmpi(paired2, 'on')
%                         disp('Removing ERSP baseline for second indep. variables (first indep. var. is unpaired)');
%                         for c = 1:nc        % ng = number of groups
%                             meanpowbase = computeerspbaseline(erspbase(c,:), opt.singletrials);
%                             ersp(c,:)   = removeerspbaseline(ersp(c,:), erspbase(c,:), meanpowbase);
%                         end;
%                     else
%                         disp('Not removing ERSP baseline (both indep. variables are unpaired');
%                     end;
%                 end;
%             end;
%         end;
%         
% %         if strcmpi(opt.statmode, 'common')
% %             % collapse the two last dimensions before computing significance
% %             % i.e. 18 subject with 4 channels -> same as 4*18 subjects
% %             % --------------------------------------------------------------
% %             disp('Using all channels for statistics...');
% %             for c = 1:nc
% %                 for g = 1:ng
% %                     ersp{c,g} = reshape( ersp{c,g}, size(ersp{c,g},1), size(ersp{c,g},2), size(ersp{c,g},3)*size(ersp{c,g},4));
% %                 end;
% %             end;
% %         end;
% 
%         % copy data to structure
%         % ----------------------
%         if ~isempty(events)
%             tmpstruct = setfield(tmpstruct, [ dtype 'events' ], events);
%         end;
%         tmpstruct = setfield(tmpstruct, [ dtype ordinate ], allfreqs);
%         tmpstruct = setfield(tmpstruct, [ dtype 'times'  ], alltimes);
%         if strcmpi(opt.singletrials, 'on')
%             tmpstruct = setfield(tmpstruct, [ dtype 'datatrials' ], ersp);
%             tmpstruct = setfield(tmpstruct, [ dtype 'subjinds' ], erspinds);
%             tmpstruct = setfield(tmpstruct, [ dtype 'times' ], alltimes);
%             if ~isempty(opt.channels)
%                  tmpstruct = setfield(tmpstruct, [ dtype 'trialinfo' ], opt.subject);
%             else tmpstruct = setfield(tmpstruct, [ dtype 'trialinfo' ], opt.component);
%             end;
%         else
%             tmpstruct = setfield(tmpstruct, [ dtype 'data' ], ersp);
%             if strcmpi(dtype, 'ersp')
%                 tmpstruct = setfield(tmpstruct, [ dtype 'base' ], erspbase);
%             end;
%         end;
%         
%         % copy results to structure
%         % -------------------------
%         fields = { [ dtype 'data'     ] [ dtype 'events' ] [ dtype ordinate ] [ dtype 'datatrials' ] ...
%                    [ dtype 'subjinds' ] [ dtype 'base'   ] [ dtype 'times'  ] [ dtype 'trialinfo'  ]  'allinds' 'setinds' };
%         for f = 1:length(fields)
%             if isfield(tmpstruct, fields{f}),
%                 tmpdata = getfield(tmpstruct, fields{f});
%                 if ~isempty(opt.channels)
%                      STUDY.changrp = setfield(STUDY.changrp, { finalinds(ind) }, fields{f}, tmpdata);
%                 else STUDY.cluster = setfield(STUDY.cluster, { finalinds(ind) }, fields{f}, tmpdata);
%                 end;
%             end;
         end;
    end;
end;

% output unit
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

% 
% if ~isempty(opt.channels)
%     structdat = STUDY.changrp;
%     erspdata  = cell(nc, ng);
%     events    = {};
%     for ind =  1:length(erspdata(:))
%         if strcmpi(opt.singletrials, 'on')
%              tmpdat = getfield(structdat(allinds(1)), [ dtype 'datatrials' ]);
%         else tmpdat = getfield(structdat(allinds(1)), [ dtype 'data' ]);
%         end;
%         if ndims(tmpdat{ind}) == 2, erspdata{ind} = zeros([ size(tmpdat{ind}) 1 length(allinds)]);
%         else                        erspdata{ind} = zeros([ size(tmpdat{ind}) length(allinds)]);
%         end;
%         for index = 1:length(allinds)
%             if strcmpi(opt.singletrials, 'on')
%                  tmpdat = getfield(structdat(allinds(index)), [ dtype 'datatrials' ]);
%             else tmpdat = getfield(structdat(allinds(index)), [ dtype 'data' ]);
%             end;
%             erspdata{ind}(:,:,:,index) = tmpdat{ind};
%             allfreqs = getfield(structdat(allinds(index)), [ dtype ordinate ]);
%             alltimes = getfield(structdat(allinds(index)), [ dtype 'times'  ]);
%             if isfield(structdat, [ dtype 'events' ])
%                 events   = getfield(structdat(allinds(index)), [ dtype 'events' ]);
%             end;
%             compinds = structdat(allinds(index)).allinds;
%             setinds  = structdat(allinds(index)).setinds;
%         end;
%         erspdata{ind} = permute(erspdata{ind}, [1 2 4 3]); % time freqs elec subjects
%     end;
%     if ~isempty(opt.subject) && strcmpi(opt.singletrials,'off')
%         erspdata = std_selsubject(erspdata, opt.subject, setinds, { STUDY.design(opt.design).cell.case }, 2); 
%     end;
% else
%     if strcmpi(opt.singletrials, 'on')
%          erspdata = getfield(STUDY.cluster(allinds(1)), [ dtype 'datatrials' ]);
%     else erspdata = getfield(STUDY.cluster(allinds(1)), [ dtype 'data' ]);
%     end;
%     allfreqs = getfield(STUDY.cluster(allinds(1)), [ dtype ordinate ]);
%     alltimes = getfield(STUDY.cluster(allinds(1)), [ dtype 'times'  ]);
%     if isfield(STUDY.cluster, [ dtype 'events' ])
%         events   = getfield(STUDY.cluster(allinds(1)), [ dtype 'events' ]);
%     end;
%     compinds = STUDY.cluster(allinds(1)).allinds;
%     setinds  = STUDY.cluster(allinds(1)).setinds;
%     if ~isempty(opt.component) && length(allinds) == 1 && strcmpi(opt.singletrials,'off')
%         erspdata = std_selcomp(STUDY, erspdata, allinds, setinds, compinds, opt.component);
%     end;
% end;

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

% call newtimef
% --------------
function [dataout tmpParams] = callnewtimef(dataSubject, xvals, yvals, pnts, tlimits, datatype, timerange, params);

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
        [tmp dataout] = newtimef([],pnts,tlimits, srate, cycles, 'verbose', 'off', params);
    end;