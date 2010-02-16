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
%  'channels'  - [cell] list of channels to import {default: all}
%  'clusters'  - [integer] list of clusters to import {[]|default: all but
%                the parent cluster (1) and any 'NotClust' clusters}
%  'timerange' - [min max] time range {default: whole measure range}
%  'subject'   - [string] select a specific subject {default:all}
%  'component' - [integer] select a specific component in a cluster
%                 {default:all}
%  'singletrials' - ['on'|'off'] load single trials spectral data (if available)
%
% Output:
%  STUDY    - updated studyset structure
%  datavals  - [cell array] erp data (the cell array size is 
%             condition x groups)
%  times    - [float array] array of time values
%  setinds  - [cell array] datasets indices
%  cinds    - [cell array] channel or component indices
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

% $Log: not supported by cvs2svn $
% Revision 1.17  2006/09/12 18:55:46  arno
% channel compatibility
%
% Revision 1.16  2006/03/29 18:23:36  scott
% help msg
%
% Revision 1.15  2006/03/14 02:27:04  scott
% help msg - filename
%
% Revision 1.14  2006/03/14 01:51:02  scott
% help msg
%
% Revision 1.13  2006/03/13 23:21:27  arno
% timerange
%
% Revision 1.12  2006/03/11 07:30:41  arno
% time input
%
% Revision 1.11  2006/03/11 07:22:22  arno
% header
%
% Revision 1.10  2006/03/11 07:18:43  arno
% header information
%
% Revision 1.9  2006/03/10 16:08:19  arno
% select time range
%
% Revision 1.8  2006/03/10 00:36:34  arno
% error msg
%
% Revision 1.7  2006/03/09 18:53:15  arno
% reading all ERPs if necessary
%
% Revision 1.6  2006/03/09 18:45:40  arno
% reading all ERP
%
% Revision 1.5  2006/03/09 18:24:36  arno
% load Matlab file now
%
% Revision 1.4  2006/03/08 20:31:25  arno
% rename func
%
% Revision 1.3  2006/03/07 22:21:26  arno
% use fullfile
%
% Revision 1.2  2006/03/07 22:09:25  arno
% fix error message
%

function [STUDY, datavals, xvals, setinds, allinds] = std_readerp(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readerp;
    return;
end
if ~isstruct(ALLEEG) % old calling format
    [STUDY, datavals] = std_readerpsub(STUDY, ALLEEG, varargin{:});
    return;
end;

STUDY = pop_erpparams(STUDY, 'default');
STUDY = pop_specparams(STUDY, 'default');
[opt moreopts] = finputcheck( varargin, { ...
    'type'          { 'string' 'cell' } { [] [] } '';
    'channels'      'cell'    []             {};
    'clusters'      'integer' []             [];
    'timerange'     'real'    []             STUDY.etc.erpparams.timerange;
    'freqrange'     'real'    []             STUDY.etc.specparams.freqrange;
    'datatype'      'string'  { 'erp' 'spec' } 'erp';
    'rmsubjmean'    'string'  { 'on' 'off' } 'off';
    'singletrials'  'string'  { 'on' 'off' } 'off';
    'component'     'integer' []             [];
    'subject'       'string'  []             '' }, ...
    'std_readerp', 'ignore');
if isstr(opt), error(opt); end;
nc = max(length(STUDY.condition),1);
ng = max(length(STUDY.group),1);
dtype = opt.datatype;

% find channel indices
% --------------------
if ~isempty(opt.channels)
     finalinds = std_chaninds(STUDY, opt.channels);
else finalinds = opt.clusters;
end;

for ind = 1:length(finalinds)

    % find indices
    % ------------
    if ~isempty(opt.channels)
        tmpstruct = STUDY.changrp(finalinds(ind));
        allinds       = tmpstruct.allinds;
        for i=1:length(allinds(:)), allinds{i} = -allinds{i}; end; % invert sign for reading
        setinds       = tmpstruct.setinds;
    else
        [ tmpstruct setinds allinds ] = std_setcomps2cell(STUDY, finalinds(ind));
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
        if isfield(tmpstruct, [ dtype 'datatrials' ]) && ~isempty(getfield(tmpstruct, [ dtype 'datatrials' ])) && eqtf
            if ~isempty(opt.channels) && strcmpi(getfield(tmpstruct, [ dtype 'trialinfo' ]), opt.subject)
                dataread = 1; 
            elseif isequal(getfield(tmpstruct, [ dtype 'trialinfo' ]), opt.component) 
                dataread = 1; 
            end;
        end;
    end;
    
    if ~dataread
        % reserve arrays
        % --------------
        alldata  = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) );
        tmpind = 1; while(isempty(setinds{tmpind})), tmpind = tmpind+1; end;
        if strcmpi(dtype, 'erp'), [ tmp xvals datatrialpresent] = std_readerpsub(ALLEEG, setinds{tmpind}(1), allinds{tmpind}(1), opt.timerange, fastif(strcmpi(opt.singletrials, 'on'),1,0));
        else                      [ tmp xvals datatrialpresent] = std_readspec(  ALLEEG, setinds{tmpind}(1), allinds{tmpind}(1), opt.freqrange, 0, fastif(strcmpi(opt.singletrials, 'on'),1,0));
        end;

        % read the data and select channels
        % ---------------------------------
        fprintf([ 'Reading ' dtype ' data...' ]);
        if strcmpi(opt.singletrials, 'on')
            if ~datatrialpresent
                fprintf('\n');
                errordlg2('No single trial data - recompute data files');
                datavals = [];
                return;
            end;
            allsubjects = { STUDY.datasetinfo(:).subject };
            for c = 1:nc
                for g = 1:ng
                    if ~isempty(opt.channels)
                        if ~isempty(opt.subject) inds = strmatch( opt.subject, allsubjects(setinds{c,g}));
                        else inds = 1:length(allinds{c,g}); end;
                    else
                        if ~isempty(opt.component) inds = find( allinds{c,g} == STUDY.cluster(finalinds(ind)).comps(opt.component));
                        else inds = 1:length(allinds{c,g}); end;
                    end;
                    if ~isempty(inds)
                        if strcmpi(dtype, 'erp') alldata{c, g} = squeeze(std_readerpsub(ALLEEG, setinds{c,g}(inds), allinds{c,g}(inds), opt.timerange, 1));
                        else                     alldata{c, g} = squeeze(std_readspec(  ALLEEG, setinds{c,g}(inds), allinds{c,g}(inds), opt.freqrange, 0, 1));
                        end;
                    end;
                end;
            end;
        else
            for c = 1:nc
                for g = 1:ng
                    if strcmpi(dtype, 'erp') alldata{c, g} = std_readerpsub( ALLEEG, setinds{c,g}(:), allinds{c,g}(:), opt.timerange)';
                    else                     alldata{c, g} = std_readspec(   ALLEEG, setinds{c,g}(:), allinds{c,g}(:), opt.freqrange)';
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
                    disp('Inverting ERP component polarities based on scalp map polarities');
                    for index = 1:length(alldata(:))
                        for comps = 1:size(alldata{index},2)
                            alldata{index}(:,comps) = alldata{index}(:,comps)*STUDY.cluster(finalinds(ind)).topopol(comps);
                        end;
                    end;
                else
                    disp('Cluster topographies absent - cannot adjust single component ERP polarities');
                end;
            end;
        end;

        % remove mean of each subject across groups and conditions
        if strcmpi(dtype, 'spec') && strcmpi(opt.rmsubjmean, 'on') && ~isempty(opt.channels) && strcmpi(opt.singletrials, 'off')
            disp('Removing mean spectrum accross subjects');
            for indtmp = 1:length(allinds{c,g}) % scan subjects
               meanspec =zeros(size( alldata{1, 1}(:,indtmp) ));
               for c = 1:nc
                    for g = 1:ng
                        meanspec = meanspec + alldata{c, g}(:,indtmp)/(nc*ng);
                    end;
               end;
               for c = 1:nc
                    for g = 1:ng
                        alldata{c, g}(:,indtmp) = alldata{c, g}(:,indtmp) - meanspec; % subtractive model
                        % alldata{c, g}(:,indtmp) = alldata{c, g}(:,indtmp)./meanspec; % divisive model
                    end;
               end;
            end;
        end;
        
        if strcmpi(opt.singletrials, 'on')
             tmpstruct = setfield( tmpstruct, [ dtype 'datatrials' ], alldata);
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
                       [ dtype 'times' ] [ dtype 'trialinfo' ] 'allinds' 'setinds' };
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

% return structure
% ----------------
allinds   = finalinds;
if ~isempty(opt.channels)
    structdat = STUDY.changrp;
    datavals = cell(nc, ng);
    for ind =  1:length(datavals(:))
        if strcmpi(opt.singletrials, 'on')
             tmpdat = getfield(structdat(allinds(1)), [ dtype 'datatrials' ]);
        else tmpdat = getfield(structdat(allinds(1)), [ dtype 'data' ]);
        end;
        datavals{ind} = zeros([ size(tmpdat{ind}) length(allinds)]);
        for index = 1:length(allinds)
            if strcmpi(opt.singletrials, 'on')
                 tmpdat = getfield(structdat(allinds(index)), [ dtype 'datatrials' ]);
            else tmpdat = getfield(structdat(allinds(index)), [ dtype 'data' ]);
            end;
            datavals{ind}(:,:,index) = tmpdat{ind};
        end;
        if strcmpi( dtype, 'spec'), xvals = getfield(structdat(allinds(index)), [ dtype 'freqs' ]);
        else                        xvals = getfield(structdat(allinds(index)), [ dtype 'times' ]);
        end;
        
        datavals{ind} = squeeze(permute(datavals{ind}, [1 3 2])); % time elec subjects
    end;
    setinds  = structdat(allinds(1)).setinds;
    if ~isempty(opt.subject) && strcmpi(opt.singletrials,'off')
        datavals = std_selsubject(datavals, opt.subject, setinds, { STUDY.datasetinfo(:).subject }, 2); 
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

% % return structure
% % ----------------
% allinds   = finalinds;
% if ~isempty(opt.channels)
%     structdat = STUDY.changrp;
%     datavals = cell(nc, ng);
%     for ind =  1:length(datavals(:))
%         if strcmpi(opt.singletrials, 'on')
%              datavals{ind} = zeros([ size(structdat(allinds(1)).datavalstrials{ind}) length(allinds)]);
%         else datavals{ind} = zeros([ size(structdat(allinds(1)).datavals{ind}) length(allinds)]);
%         end;
%         for index = 1:length(allinds)
%             if strcmpi(opt.singletrials, 'on')
%                  datavals{ind}(:,:,index) = structdat(allinds(index)).datavalstrials{ind};
%             else datavals{ind}(:,:,index) = structdat(allinds(index)).datavals{ind};
%             end;
%             xvals                 = structdat(allinds(index)).erptimes;
%             compinds                 = structdat(allinds(index)).allinds;
%             setinds                  = structdat(allinds(index)).setinds;
%         end;
%         datavals{ind} = squeeze(permute(datavals{ind}, [1 3 2])); % time elec subjects
%         datavals{ind} = squeeze(permute(datavals{ind}, [1 3 2])); % time elec subjects
%     end;
%     if ~isempty(opt.subject) && strcmpi(opt.singletrials,'off')
%         datavals = std_selsubject(datavals, opt.subject, setinds, { STUDY.datasetinfo(:).subject }, 2); 
%     end;
% else
%     if strcmpi(opt.singletrials, 'on')
%          datavals = STUDY.cluster(allinds(1)).datavalstrials;
%     else datavals = STUDY.cluster(allinds(1)).datavals;
%     end;
%     xvals = STUDY.cluster(allinds(1)).erptimes;
%     compinds = STUDY.cluster(allinds(1)).allinds;
%     setinds  = STUDY.cluster(allinds(1)).setinds;
%     if ~isempty(opt.component) && length(allinds) == 1 && strcmpi(opt.singletrials,'off')
%         datavals = std_selcomp(STUDY, datavals, allinds, setinds, compinds, opt.component);
%     end;
% end;

% std_readerp() - returns the ERP for an ICA component in an epoched dataset.
%                 The ERPs of the dataset ICA components are assumed to have 
%                 been saved in a Matlab file, [dataset_name].icaerp, in the
%                 same directory as the dataset file. If this file doesn't exist, 
%                 use std_erp() to create it, else use a pre-clustering function 
%                 that calls it: pop_preclust() or std_preclust()  
% Usage:    
%   >> [erp, times] = std_readerp(ALLEEG, setindx, component, timewindow);  
%
% Inputs:
%   ALLEEG     - an EEGLAB dataset vector (else one EEG dataset). 
%                ALLEEG must contain the dataset of interest (see 'setindx').
%   setindx    - [integer] index of the EEG dataset in the ALLEEG structure 
%                for which to read the component ERP.
%   component  - [integer] index of the component in the selected EEG dataset 
%                for which to return the ERP. 
%   timewindow - [min max] ERP time (latency) window, in ms. Must be in
%                the dataset epoch latency range.
% Outputs:
%   erp        - ERP for the requested ICA component in the selected dataset; 
%                the average of the ICA activations in all the dataset epochs.
%   times      - vector of ERP time points (latencies) in ms.
%
%  See also  std_erp(), pop_preclust(), std_preclust()          
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, February, 2005

function [X, t, singletrialdatapresent] = std_readerpsub(ALLEEG, abset, comp, timerange, singletrial)

if nargin < 4
    timerange = [];
end;
if nargin < 5
    singletrial = 0;
end;

X = [];
if length(abset) < length(comp)
    abset = ones(1,length(comp))*abset;
end;

% convert components or channel indices
% -------------------------------------
if iscell(comp)
    % find channel indices list
    % -------------------------
    chanind  = [];
    chanlabs = lower({ ALLEEG(abset(1)).chanlocs.labels });
    for index = 1:length(comp)
        tmp = strmatch(lower(comp{index}), chanlabs, 'exact');
        if isempty(tmp)
            error([ 'Channel ''' comp{index} ''' not found in dataset ' int2str(abset)]);
        else    
            chanind = [ chanind tmp ];
        end;
    end;
    filename = fullfile( ALLEEG(abset(1)).filepath,[ ALLEEG(abset(1)).filename(1:end-3) 'daterp']);
    prefix = 'chan';
    inds   = chanind;
elseif comp(1) < 0
    filename = fullfile( ALLEEG(abset(1)).filepath,[ ALLEEG(abset(1)).filename(1:end-3) 'daterp']);
    prefix = 'chan';
    inds   = -comp;
else
    filename = fullfile( ALLEEG(abset(1)).filepath,[ ALLEEG(abset(1)).filename(1:end-3) 'icaerp']);
    prefix = 'comp';
    inds   = comp;
end;

singletrialdatapresent = 1;
for k = 1:length(abset)

    if strcmpi(prefix, 'chan')
         filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'daterp']);
    else filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'icaerp']);
    end;

    try,
        erpstruct = load( '-mat', filename, [ prefix int2str(inds(k)) ], 'times' );
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end;
    
    tmpdat    = getfield(erpstruct, [ prefix int2str(inds(k)) ]);
    if singletrial == 0,
        if size(tmpdat,2) > 1 && size(tmpdat,1) > 1, tmpdat = mean(tmpdat,2); end;
        if k == 1
            if size(tmpdat,1) == 1, tmpdat = tmpdat'; end;
            X = zeros([ length(comp) size(tmpdat) ]);
        end;
        X(k,:,:) = tmpdat;
    else
        if size(tmpdat,1) == 1 || size(tmpdat,2) == 1
            singletrialdatapresent = 0;
        end;
        if k == 1
            if size(tmpdat,1) == 1, tmpdat = tmpdat'; end;
            X = zeros([ 1 size(tmpdat) ]);
            X(1,:,:) = tmpdat;
        else
            X(1,:,end+1:end+size(tmpdat,2)) = tmpdat;
        end;
    end;
    t = getfield(erpstruct, 'times');
end;

% select time range of interest
% -----------------------------
if ~isempty(timerange)
    maxind = max(find(t <= timerange(end)));
    minind = min(find(t >= timerange(1)));
else
    maxind = length(t);
    minind = 1;
end;
X = X(:,minind:maxind,:);
t = t(minind:maxind)';
