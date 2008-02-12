% std_readdata() - load one or more requested measures 
%                  ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map']
%                  for all components of a specified cluster.  
%                  Called by cluster plotting functions 
%                  std_envtopo(), std_erpplot(), std_erspplot(), ...
% Usage:
%         >> [STUDY, clustinfo, finalinds] = std_readdata(STUDY, ALLEEG); 
%                                              % return all measures
%         >> [STUDY, clustinfo, finalinds] = std_readdata(STUDY, ALLEEG, ...
%                                              cluster, infotype,varargin);
% Inputs:
%       STUDY - studyset structure containing some or all files in ALLEEG
%      ALLEEG - vector of loaded EEG datasets
%
% Optional inputs:
%  'infotype'  - ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map'|'data'|'event'] 
%                type of stored data or cluster information to read. May 
%                also be a cell array of these types, for example: { 'erp' 
%                'map' 'dipole' }. {default: 'erp'}
%  'channels'  - [cell] list of channels to import {default: all}
%  'clusters'  - [integer] list of clusters to import {[]|default: all but
%                the parent cluster (1) and any 'NotClust' clusters}
%  'freqrange' - [min max] frequency range {default: whole measure range}
%  'timerange' - [min max] time range {default: whole measure epoch}
%  'statmode'  - ['individual'|'trials'] statistical mode for ERSP (also 
%                specify what type of data to import -- mean (of individual 
%                subjects), or trials. This functionality is still unstable
%                for 'trials' { default: 'individual'}
%  'rmsubjmean' - ['on'|'off'] remove mean subject spectrum from every component
%                 or channel spectrum, making them easier to compare
%                 { default: 'off' }
%  'subbaseline' - ['on'|'off'] remove all condition and spectrum baselines
%                  for ERSP data { default: 'on' }
%
% Optional inputs for events:
%  'type','timewin','fieldname' - optional parameters for the function 
%                  eeg_geteventepoch may also be used to select specific
%                  events.
%
% Optional inputs for raw data:
%  'rmclust'  - [integer] list of artifact cluster indices to subtract 
%               when reading raw data.
%
% Output:
%    STUDY     - (possibly) updated STUDY structure
%    clustinfo - structure of specified cluster information. 
%                This is the same as STUDY.cluster(cluster_number)     
%    Fields:
%     clustinfo.erpdata    % (ncomps, ntimes) array of component ERPs
%     clustinfo.erptimes   % vector of ERP epoch latencies (ms)
%
%     clustinfo.specdata   % (ncomps, nfreqs) array of component spectra
%     clustinfo.specfreqs  % vector of spectral frequencies (Hz)
%
%     clustinfo.erspdata   % (ncomps,ntimes,nfreqs) array of component ERSPs
%      clustinfo.ersptimes % vector of ERSP latencies (ms)
%      clustinfo.erspfreqs % vector of ERSP frequencie (Hz)
%
%     clustinfo.itcdata    % (ncomps,ntimes,nfreqs) array of component ITCs
%      clustinfo.itctimes  % vector of ITC latencies (ms)
%      clustinfo.itc_freqs % vector of ITC frequencies (Hz)
%
%     clustinfo.topo       % (ncomps,65,65) array of component scalp map grids
%       clustinfo.topox    % abscissa values for columns of the scalp maps
%       clustinfo.topoy    % ordinate values for rows of the scalp maps
%
%     clustinfo.data           % raw data arrays
%       clustinfo.datatimes    % time point for raw data
%       clustinfo.datasortvals % sorting values (one per trial)
%       clustinfo.datacontinds % dataset indices (contacted for all trials)
%
%     clustinfo.dipole     % array of component dipole information structs
%                          % with the same format as EEG.dipfit.model
%
%   finalinds - either the cluster(s) or channel(s) indices selected.
%
% Example:
%         % To obtain the ERPs for all Cluster-3 components from a STUDY
%         %
%         [STUDY clustinfo] = std_readdata(STUDY, ALLEEG, 'clusters',3, 'infotype','erp');
%         figure; plot(clustinfo.erptimes, mean(clustinfo.erpdata,2));
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
% Revision 1.41  2008/02/12 19:07:20  nima
% _
%
% Revision 1.40  2008/02/08 20:21:42  arno
% *** empty log message ***
%
% Revision 1.39  2008/01/10 23:58:14  arno
% defining event sorting fields explicitelly
%
% Revision 1.38  2008/01/10 20:00:10  arno
% update help message
%
% Revision 1.37  2008/01/10 19:56:18  arno
% adding even reading
%
% Revision 1.36  2008/01/08 23:05:41  arno
% added data type to read raw data
%
% Revision 1.35  2007/09/11 10:51:33  arno
% fix indices problems preventing crashes
%
% Revision 1.34  2007/08/23 01:59:36  arno
% remove debugging comment
%
% Revision 1.33  2007/08/23 01:58:38  arno
% fix nc and ng
%
% Revision 1.32  2007/08/13 18:01:36  arno
% remove all reference to chaninds
%
% Revision 1.31  2007/08/08 00:31:17  nima
% the help was incorrect, now we have to explicitly specify the key values, so 'clusters',3, 'infotype','erp'
%
% Revision 1.30  2007/08/06 23:09:55  scott
% help msg clarification - please check defaults, Arno.
% In particular, why is the default for 'infotype' = 'erp' instead of
% (new option) 'all' ?
% Is the STUDY.cluster structure returned in the 2nd arg ?
% Can this function also return the IC and subject numbers for the cluster ICs??
%
% Revision 1.29  2007/08/06 22:23:20  arno
% remove unused code
%
% Revision 1.28  2007/08/06 22:01:53  arno
% implement map and dipoles
%
% Revision 1.27  2007/08/06 21:49:40  arno
% reading dipoles
% /
%
% Revision 1.26  2007/08/06 21:18:54  arno
% fix header
%
% Revision 1.25  2007/07/27 22:22:48  toby
% nothing
%
% Revision 1.24  2007/06/08 20:10:41  toby
% Help updated to state optional inputs and the default values. 'subbaseline' and 'rmsubmean' defaults reset to 'off' because they were causing crashes with some studies. It was not clear how they were intended to function, so I did not fix them.
%
% Revision 1.23  2007/06/02 03:12:21  toby
% trying to fix a NaN indexing problem
%
% Revision 1.22  2007/05/11 02:09:23  toby
% bug when ersp/itc of nonexistent component/channel queried fixed (poorly)
%
% Revision 1.21  2007/04/05 21:33:31  arno
% same
%
% Revision 1.20  2007/04/05 21:31:58  arno
% same
%
% Revision 1.19  2007/04/05 21:30:54  arno
% conversion to single, Matlab 6.5 compatibility
%
% Revision 1.18  2007/03/17 21:17:59  arno
% fix if condition string is a substring of another condition string
%
% Revision 1.17  2007/01/26 18:06:15  arno
% minor bug
%
% Revision 1.16  2006/11/23 01:12:46  arno
% implement mean spectrum subtraction
%
% Revision 1.15  2006/11/15 21:52:14  arno
% reading baseline ERSP if spectrum absent
%
% Revision 1.14  2006/11/09 22:38:08  arno
% fix spectrum loading
%
% Revision 1.13  2006/11/09 20:50:34  arno
% final indices now returned
%
% Revision 1.12  2006/11/09 01:32:28  arno
% new implementation
%
% Revision 1.11  2006/11/09 00:35:59  arno
% text if condgrp exist
%
% Revision 1.10  2006/11/09 00:27:01  arno
% same
%
% Revision 1.9  2006/11/09 00:20:38  arno
% same
%
% Revision 1.8  2006/11/09 00:16:44  arno
% field compsind
%
% Revision 1.7  2006/11/08 23:46:48  arno
% save condgrp structure
%
% Revision 1.6  2006/11/02 22:05:16  arno
% same
%
% Revision 1.5  2006/11/02 22:03:35  arno
% add subject statmode
%
% Revision 1.4  2006/09/20 12:27:14  arno
% fix reading of spectral data
%

function [STUDY, clustinfo, finalinds] = std_readdata(STUDY, ALLEEG, varargin);

if nargin < 2
    help std_readdata;
    return;
end
if nargin < 3
    channel = {};
end

[opt moreopts] = finputcheck( varargin, { ...
    'type'       { 'string' 'cell' } { [] [] } '';
    'timewin'    'real'    []        [-Inf Inf];
    'fieldname'  'string'  []        'latency';
    'condition'  'cell'    []       {};
    'channels'   'cell'    []       {};
    'clusters'   'integer' []       [];
    'rmclust'    'integer' []       [];
    'freqrange'  'real'    []       [];
    'timerange'  'real'    []       [];
    'statmode'   'string'  { 'subjects' 'individual' 'common' 'trials' }       'individual';
    'subbaseline' 'string'  { 'on' 'off' }       'off';
    'rmsubjmean'  'string'  { 'on' 'off' }       'off';
    'infotype'   'string'  { 'erp' 'spec' 'ersp' 'itc' 'map' 'topo' 'dipole' 'scalp' 'data' 'event' } 'erp' }, ...
                              'std_readdata', 'ignore');
if isstr(opt), error(opt); end;

if strcmpi(opt.infotype, 'erp'),
    STUDY = pop_erpparams(STUDY, 'default');
    if isempty(opt.timerange), opt.timerange = STUDY.etc.erpparams.timerange; end;
elseif strcmpi(opt.infotype, 'spec'),
    STUDY = pop_specparams(STUDY, 'default');
    if isempty(opt.freqrange), opt.freqrange = STUDY.etc.specparams.freqrange; end;
elseif strcmpi(opt.infotype, 'ersp') | strcmpi(opt.infotype, 'itc')
    STUDY = pop_erspparams(STUDY, 'default');
    if isempty(opt.freqrange), opt.freqrange   = STUDY.etc.erspparams.freqrange; end;
    if isempty(opt.timerange), opt.timerange   = STUDY.etc.erspparams.timerange; end;
    if strcmpi(opt.statmode, 'individual') |  strcmpi(opt.statmode, 'subjects'), opt.statmode    = STUDY.etc.erspparams.statmode;end;
    if strcmpi(opt.subbaseline, 'on'),      opt.subbaseline = STUDY.etc.erspparams.subbaseline; end;
end;

nc = max(length(STUDY.condition),1);
ng = max(length(STUDY.group),1);
zero = single(0);
tmpver = version;
if tmpver(1) == '6' | tmpver(1) == '5'
    zero = double(zero);
end;

% find channel indices
% --------------------
if ~isempty(opt.channels)
    finalinds = std_chaninds(STUDY, opt.channels);
else
    finalinds = opt.clusters;
end;
    
% read topography with another function
% -------------------------------------
if strcmpi(opt.infotype, 'map') | strcmpi(opt.infotype, 'scalp') | strcmpi(opt.infotype, 'topo')
    [STUDY tmpclust] = std_readtopoclust(STUDY, ALLEEG, opt.clusters);
    clustinfo = [];
    for index = 1:length(tmpclust)
        if index == 1, clustinfo        = tmpclust{index};
        else           clustinfo(index) = tmpclust{index};
        end;
    end;
    return;
end;

% read indices for removing component clusters from data
% ------------------------------------------------------
if ~isempty(opt.rmclust)
    for ind = 1:length(opt.rmclust)
        [tmp setindsrm{ind} allindsrm{ind}] = getsetinds(STUDY, opt.rmclust(ind));
    end;
end;

for ind = 1:length(finalinds)
    
    % find indices
    % ------------
    if ~isempty(opt.channels)
         tmpstruct = STUDY.changrp(finalinds(ind));
         alldatasets = 1:length(STUDY.datasetinfo);
         %allchanorcomp = -tmpstruct.chaninds;
         allinds       = tmpstruct.allinds;
         for i=1:length(allinds(:)), allinds{i} = -allinds{i}; end; % invert sign
         setinds       = tmpstruct.setinds;
    else 
        [ tmpstruct setinds allinds ] = getsetinds(STUDY, finalinds(ind));
    end;
                    
    dataread = 0;
    switch opt.infotype
        case 'event',
            % check if data is already here
            % -----------------------------
            if isfield(tmpstruct, 'datasortvals')
                if ~isempty(tmpstruct.datasortvals)
                    dataread = 1;
                end;
            end;
            
            if ~dataread
                % reserve arrays
                % --------------
                datasortvals   = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) );
                datacontinds    = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) );
                for c = 1:nc
                    for g = 1:ng
                        datasortvals{c, g} = repmat(zero, [1, sum([ ALLEEG(setinds{c,g}).trials ] )]);
                        datacontinds{c, g} = repmat(zero, [1, sum([ ALLEEG(setinds{c,g}).trials ] )]);
                    end;
                end;
            
                % read the data and select channels
                % ---------------------------------
                if ~isempty(opt.type)
                    fprintf('Reading events:');
                    for c = 1:nc
                        for g = 1:ng
                            counttrial = 1;
                            for indtmp = 1:length(allinds{c,g})
                                tmpdata  = eeg_getepochevent(ALLEEG(setinds{c,g}(indtmp)), ...
                                                             'type', opt.type, 'timewin', opt.timewin, 'fieldname', opt.fieldname);
                                datasortvals{c, g}(:,counttrial:counttrial+ALLEEG(setinds{c,g}(indtmp)).trials-1) = squeeze(tmpdata);
                                datacontinds{c, g}(:,counttrial:counttrial+ALLEEG(setinds{c,g}(indtmp)).trials-1) = setinds{c,g}(indtmp);
                                counttrial = counttrial+ALLEEG(setinds{c,g}(indtmp)).trials;
                                fprintf('.');
                            end;
                        end;
                    end;
                end;
                fprintf('\n');
                tmpstruct.datasortvals = datasortvals;
                tmpstruct.datacontinds = datacontinds;
            end;
     
        case 'data',
            % check if data is already here
            % -----------------------------
            if isfield(tmpstruct, 'data')
                if ~isempty(tmpstruct.data)
                    dataread = 1;
                end;
            end;
            
            if ~dataread
                % reserve arrays
                % --------------
                alldata   = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) );
                for c = 1:nc
                    for g = 1:ng
                        alldata{c, g} = repmat(zero, [ALLEEG(1).pnts, sum([ ALLEEG(setinds{c,g}).trials ] )]);
                    end;
                end;
            
                % read the data and select channels
                % ---------------------------------
                fprintf('Reading data/activity:');
                for c = 1:nc
                    for g = 1:ng
                        counttrial = 1;
                        for indtmp = 1:length(allinds{c,g})
                            settmpind = STUDY.datasetinfo(setinds{c,g}(indtmp)).index;
                            if ~isempty(opt.channels)
                                tmpdata = eeg_getdatact(ALLEEG(settmpind),  'channel', -allinds{c,g}(indtmp), 'verbose', 'off');
                            else
                                tmpdata = eeg_getdatact(ALLEEG(settmpind), 'component', allinds{c,g}(indtmp), 'verbose', 'off');
                            end;
                            
                            % remove artifactual components
                            % -----------------------------
                            if ~isempty(opt.rmclust)                            
                                for rmi = 1:length(opt.rmclust)
                                    findind   = find(settmpind == STUDY.cluster(opt.rmclust(rmi)).setinds{c,g});
                                    rmcomps   = [ rmcomps allindsrm{c,g}(findind) ];
                                end;
                                ws = ALLEEG(settmpind).icaweights * ALLEEG(settmpind).icasphere;
                                w = ALLEEG(settmpind).icawinv(-allinds{c,g}(indtmp),rmcomps)*ws(rmcomps,:);
                                %tmpdata = tmpdata - 
                            end;
                                
                            alldata{c, g}(:,counttrial:counttrial+ALLEEG(setinds{c,g}(indtmp)).trials-1) = squeeze(tmpdata);
                            counttrial = counttrial+ALLEEG(setinds{c,g}(indtmp)).trials;
                            fprintf('.');
                        end;
                    end;
                end;
                fprintf('\n');
                tmpstruct.datatimes = ALLEEG(1).times;
                tmpstruct.data      = alldata;
            end;
        
        case 'erp', 
            % check if data is already here
            % -----------------------------
            if isfield(tmpstruct, 'erpdata')
                if isequal( STUDY.etc.erpparams.timerange, opt.timerange) & ~isempty(tmpstruct.erpdata)
                    dataread = 1;
                end;
            end;

            if ~dataread
                % reserve arrays
                % --------------
                allerp   = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) );
                tmpind = 1; while(isempty(setinds{tmpind})), tmpind = tmpind+1; end;
                [ tmp alltimes ] = std_readerp( ALLEEG, setinds{tmpind}(1), allinds{tmpind}(1), opt.timerange);
                for c = 1:nc
                    for g = 1:ng
                        allerp{c, g} = repmat(zero, [length(alltimes), length(allinds{c,g})]);
                    end;
                end;

                % read the data and select channels
                % ---------------------------------
                fprintf('Reading ERP data:');
                for c = 1:nc
                    for g = 1:ng
                        for indtmp = 1:length(allinds{c,g})
                            [ tmperp alltimes ] = std_readerp( ALLEEG, setinds{c,g}(indtmp), allinds{c,g}(indtmp), opt.timerange);
                            allerp{c, g}(:,indtmp) = tmperp(:);
                            fprintf('.');
                        end;
                    end;
                end;
                fprintf('\n');
                tmpstruct.erptimes = alltimes;
                tmpstruct.erpdata  = allerp;
            end;

        case 'spec', 
            % check if data is already here
            % -----------------------------
            if isfield(tmpstruct, 'specdata')
                if isequal( STUDY.etc.specparams.freqrange, opt.freqrange) & ~isempty(tmpstruct.specdata)
                    dataread = 1;
                end;
            end;

            if ~dataread
                % reserve arrays
                % --------------
                allspec  = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) ); 
                filetype = 'spec';
                tmpind = 1; while(isempty(setinds{tmpind})), tmpind = tmpind+1; end;
                try,
                    [ tmp allfreqs ] = std_readspec( ALLEEG, setinds{tmpind}(1), allinds{tmpind}(1), opt.freqrange);
                catch 
                    filetype = 'ersp';
                    disp('Cannot find spectral file, trying ERSP baseline file instead');
                    [ tmpersp allfreqs alltimes tmpparams tmpspec] = std_readersp( ALLEEG, setinds{tmpind}(1), allinds{tmpind}(1), [], opt.freqrange);
                end;
                for c = 1:nc
                    for g = 1:ng
                        allspec{c, g} = repmat(zero, [length(allfreqs), length(allinds{c,g}) ]);
                    end;
                end;

                % read the data and select channels
                % ---------------------------------
                fprintf('Reading Spectrum data:');
                for c = 1:nc
                    for g = 1:ng
                        for indtmp = 1:length(allinds{c,g})
                            if strcmpi(filetype, 'spec')
                                [ tmpspec allfreqs ] = std_readspec( ALLEEG, setinds{c,g}(indtmp), allinds{c,g}(indtmp), opt.freqrange, strcmpi(opt.rmsubjmean, 'on'));
                                allspec{c, g}(:,indtmp) = tmpspec(:);
                            else
                                [ tmpersp allfreqs alltimes tmpparams tmpspec] = std_readersp( ALLEEG, setinds{c,g}(indtmp), allinds{c,g}(indtmp), [], opt.freqrange);
                                allspec{c, g}(:,indtmp) = 10*log(tmpspec(:));
                            end;
                            fprintf('.');
                        end;
                    end;
                end;
                fprintf('\n');
                tmpstruct.specfreqs = allfreqs;
                tmpstruct.specdata  = allspec;
            end;
            
        case { 'ersp' 'itc' }, 
            % check if data is already here
            % -----------------------------
            if strcmpi(opt.infotype, 'ersp')
                if isfield(tmpstruct, 'erspdata')
                    if isequal( STUDY.etc.erspparams.timerange, opt.timerange) & ...
                       isequal( STUDY.etc.erspparams.freqrange, opt.freqrange) & ~isempty(tmpstruct.erspdata)
                        dataread = 1;
                    end;
                end;
            else
                if isfield(tmpstruct, 'itcdata')
                    if isequal( STUDY.etc.erspparams.timerange, opt.timerange) & ...
                       isequal( STUDY.etc.erspparams.freqrange, opt.freqrange) & ~isempty(tmpstruct.itcdata)
                        dataread = 1;
                    end;
                end;
            end;
            
            if ~dataread
                % find total nb of trials
                % -----------------------
                if strcmpi(opt.statmode, 'trials')
                    tottrials = cell( length(STUDY.condition), length(STUDY.group) );
                    for index = 1:length( STUDY.datasetinfo )
                        condind = strmatch( STUDY.datasetinfo(index).condition, STUDY.condition );
                        grpind  = strmatch( STUDY.datasetinfo(index).group    , STUDY.group     );
                        if isempty(tottrials{condind, grpind}), tottrials{condind, grpind} = ALLEEG(index).trials;
                        else       tottrials{condind, grpind} = tottrials{condind, grpind} + ALLEEG(index).trials;
                        end;
                    end;
                end;
                
                % reserve arrays
                % --------------
                ersp     = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) );
                erspbase = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) );
                erspinds = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) );
                tmpind = 1; while(isempty(setinds{tmpind})), tmpind = tmpind+1; end;
                [ tmp allfreqs alltimes ] = std_readersp( ALLEEG, setinds{tmpind}(1), Inf*allinds{tmpind}(1), opt.timerange, opt.freqrange);
                for c = 1:nc
                    for g = 1:ng
                        ersp{c, g}     = repmat(zero, [length(alltimes), length(allfreqs), length(allinds{c,g}) ]);
                        erspbase{c, g} = repmat(zero, [               1, length(allfreqs), length(allinds{c,g}) ]);
                        if strcmpi(opt.statmode, 'trials')
                            ersp{ c, g} = repmat(zero, [length(alltimes), length(allfreqs), tottrials{c, g} ]);
                            count{c, g} = 1;
                        end;
                    end;
                end;
            
                % read the data and select channels
                % ---------------------------------
                fprintf('Reading all %s data:', upper(opt.infotype));
                for c = 1:nc
                    for g = 1:ng
                        for indtmp = 1:length(allinds{c,g})
                            try     % this 'try' is a poor solution the problem of attempting 
                                    % to read specific channel/component data that doesn't exist
                                    % called below by: allinds{c,g}(indtmp)
                            if strcmpi(opt.statmode, 'trials')
                                [ tmpersp allfreqs alltimes tmpparams] = std_readtimef( ALLEEG, setinds{c,g}(indtmp), allinds{c,g}(indtmp), ...
                                                                                  opt.timerange, opt.freqrange);
                                indices = [count{c, g}:count{c, g}+size(tmpersp,3)-1];
                                ersp{    c, g}(:,:,indices) = permute(tmpersp, [2 1 3]);
                                erspinds{c, g}(1:2,indtmp) = [ count{c, g} count{c, g}+size(tmpersp,3)-1 ];
                                count{c, g} = count{c, g}+size(tmpersp,3);
                                if size(tmpersp,3) ~= ALLEEG(STUDY.datasetinfo(index).index).trials
                                    error( sprintf('Wrong number of trials in datafile for dataset %d\n', STUDY.datasetinfo(index).index));
                                end;
                            elseif strcmpi(opt.infotype, 'itc')
                                [ tmpersp allfreqs alltimes tmpparams] = std_readitc( ALLEEG, setinds{c,g}(indtmp), allinds{c,g}(indtmp), ...
                                                                                  opt.timerange, opt.freqrange);
                                ersp{c, g}(:,:,indtmp)     = abs(permute(tmpersp    , [2 1]));
                            else
                                [ tmpersp allfreqs alltimes tmpparams tmperspbase] = std_readersp( ALLEEG, setinds{c,g}(indtmp), allinds{c,g}(indtmp), ...
                                                                                  opt.timerange, opt.freqrange);
                                ersp{c, g}(    :,:,indtmp) = permute(tmpersp    , [2 1]);
                                erspbase{c, g}(:,:,indtmp) = 10*log(permute(tmperspbase, [2 1]));
                            end
                            catch
                            end
                            fprintf('.');
                        end;
                    end;
                end;
                fprintf('\n');

                % compute ERSP or ITC if trial mode
                % (since only the timef have been loaded)
                % ---------------------------------------
                if strcmpi(opt.statmode, 'trials')
                    for c = 1:nc
                        for g = 1:ng
                            if strcmpi(opt.infotype, 'itc')
                                ersp{c,g} = ersp{c,g}./abs(ersp{c,g});
                            else
                                ersp{c,g} = 20*log10(abs(ersp{c,g}));
                                % remove baseline (each trial baseline is removed => could
                                % also be the average across all data trials)
                                [tmp indl] = min( abs(alltimes-0) );
                                erspbase{c,g} = mean(ersp{c,g}(1:indl,:,:,:));
                                ersp{c,g} = ersp{c,g} - repmat(erspbase{c,g}, [size(ersp{c,g},1) 1 1 1]);
                            end;
                        end;
                    end;
                end;
                
                % compute average baseline across groups and conditions
                % -----------------------------------------------------
                if strcmpi(opt.subbaseline, 'on') & strcmpi(opt.infotype, 'ersp')
                    disp('Recomputing baseline...');
                    for g = 1:ng        % ng = number of groups
                        for c = 1:nc    % nc = number of components
                            if strcmpi(opt.statmode, 'trials')
                                if c == 1, meanpowbase = abs(mean(erspbase{c,g}/nc,3));
                                else       meanpowbase = meanpowbase + abs(mean(erspbase{c,g}/nc,3));
                                end;
                            else
                                if c == 1, meanpowbase = abs(erspbase{c,g}/nc);
                                else       meanpowbase = meanpowbase + abs(erspbase{c,g}/nc);
                                end;
                            end;
                        end;
 
                        % subtract average baseline
                        % -------------------------
                        for c = 1:nc
                            if strcmpi(opt.statmode, 'trials'), tmpmeanpowbase = repmat(meanpowbase, [length(alltimes) 1 tottrials{c,g}]);
                            else                                tmpmeanpowbase = repmat(meanpowbase, [length(alltimes) 1 1]);
                            end;
                            ersp{c,g} = ersp{c,g} - repmat(abs(erspbase{c,g}), [length(alltimes) 1 1 1]) + tmpmeanpowbase;
                        end;
                        clear meanpowbase;
                    end;
                end;

                if strcmpi(opt.statmode, 'common')
                    % collapse the two last dimensions before computing significance
                    % i.e. 18 subject with 4 channels -> same as 4*18 subjects
                    % --------------------------------------------------------------
                    disp('Using all channels for statistics...');
                    for c = 1:nc
                        for g = 1:ng
                            ersp{c,g} = reshape( ersp{c,g}, size(ersp{c,g},1), size(ersp{c,g},2), size(ersp{c,g},3)*size(ersp{c,g},4));
                        end;
                    end;
                end;
                
                % copy data to structure
                % ----------------------
                if strcmpi(opt.infotype, 'ersp')
                    tmpstruct.erspfreqs = allfreqs;
                    tmpstruct.ersptimes = alltimes;
                    tmpstruct.erspdata  = ersp;
                    tmpstruct.erspbase  = erspbase;
                    if strcmpi(opt.statmode, 'trials')
                        tmpstruct.erspsubjinds  = erspinds;
                    end;
                else
                    tmpstruct.itcfreqs = allfreqs;
                    tmpstruct.itctimes = alltimes;
                    tmpstruct.itcdata  = ersp;
                    if strcmpi(opt.statmode, 'trials')
                        tmpstruct.itcsubjinds  = erspinds;
                    end;
                end;
            end;
        case 'dipole',                        
         fprintf('Reading dipole data...\n');
         alldips = {};
         for c = 1:nc
             for g = 1:ng
                 for indtmp = 1:length(allinds{c,g})
                     alldips{c, g}(indtmp) = ALLEEG(setinds{c,g}(indtmp)).dipfit.model(allinds{c,g}(indtmp));
                 end;
             end;
         end;
         tmpstruct.dipoles = alldips;

        case { 'map' 'scalp' 'topo' }
            % this is currenlty being done by the function std_readtopoclust
            % at the beginning of this function
        otherwise, error('Unrecognized ''infotype'' entry');
    end; % end switch
    
    % copy results to structure
    % -------------------------
    fieldnames = { 'erpdata' 'erptimes' 'specdata' 'specfreqs' 'erspdata' 'erspbase' 'erspfreqs' 'ersptimes' ...
                   'itcfreqs' 'itctimes' 'itcdata' 'erspsubjinds' 'itcsubjinds' 'allinds' 'setinds' 'dipoles' ...
                   'data' 'datatimes' 'datasortvals' 'datacontinds' };
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

% return structure
% ----------------
if ~isempty(opt.channels)
     clustinfo = STUDY.changrp(finalinds);
else clustinfo = STUDY.cluster(finalinds);
end;

% get set and indices for components
% ----------------------------------
function [ tmpstruct setinds allinds ] = getsetinds(STUDY, ind)

    tmpstruct = STUDY.cluster(ind);
    alldatasets = tmpstruct.sets; 
    allchanorcomp = repmat(tmpstruct.comps, [1 size(alldatasets,1)]);
         
    alldatasets   = alldatasets(:)';
    allchanorcomp = allchanorcomp(:)';
    
    % get indices for all groups and conditions
    % -----------------------------------------
    nc = length(STUDY.condition);
    ng = length(STUDY.group);
    allinds = cell( nc, ng );
    setinds = cell( nc, ng );
    for indtmp = 1:length(alldatasets)
        if ~isnan(alldatasets(indtmp))
            index = alldatasets(indtmp);
            condind = strmatch( STUDY.datasetinfo(index).condition, STUDY.condition, 'exact'); if isempty(condind), condind = 1; end;
            grpind  = strmatch( STUDY.datasetinfo(index).group    , STUDY.group    , 'exact'); if isempty(grpind) , grpind  = 1; end;
            indcellarray = length(allinds{condind, grpind})+1;
        end
        % load data
        % ---------
        tmpind = allchanorcomp(indtmp); 
        if ~isnan(tmpind)
            allinds{ condind, grpind}(indcellarray) = tmpind;                    
            setinds{ condind, grpind}(indcellarray) = index;                    
        end;
    end;
    tmpstruct.allinds = allinds;
    tmpstruct.setinds = setinds;

