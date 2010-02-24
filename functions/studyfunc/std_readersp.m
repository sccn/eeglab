% std_readersp() - load spectrum measures for data channels or 
%                  for all components of a specified cluster.
%                  Called by plotting functions
%                  std_envtopo(), std_erpplot(), std_erspplot(), ...
% Usage:
%         >> [STUDY, specdata, allfreqs, setinds, cinds] = ...
%                   std_readersp(STUDY, ALLEEG, varargin);
% Inputs:
%       STUDY - studyset structure containing some or all files in ALLEEG
%      ALLEEG - vector of loaded EEG datasets
%
% Optional inputs:
%  'channels'  - [cell] list of channels to import {default: all}
%  'clusters'  - [integer] list of clusters to import {[]|default: all but
%                the parent cluster (1) and any 'NotClust' clusters}
%  'freqrange' - [min max] frequency range {default: whole measure range}
%  'subject'    - [string] select a specific subject {default:all}
%  'component'  - [integer] select a specific component in a cluster
%                 {default:all}
%  'singletrials' - ['on'|'off'] load single trials data (if available)
%
% Output:
%  STUDY    - updated studyset structure
%  erspdata - [cell array] ERSP data (the cell array size is 
%             condition x groups)
%  freqs    - [float array] array of frequencies
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
% Revision 1.36  2010/02/16 08:43:21  arno
% New single-trial reading/writing
%
% Revision 1.35  2009/05/12 17:53:00  arno
% fixing dimension problem
%
% Revision 1.34  2009/05/12 17:47:20  arno
% fix baseline orientation
%
% Revision 1.33  2009/05/11 22:24:06  arno
% fix orientation problem for ERSP baseline
%
% Revision 1.32  2007/10/25 21:18:28  nima
% output parameters updated.
%
% Revision 1.31  2007/08/06 19:32:18  arno
% fix last changes
%
% Revision 1.30  2007/08/06 19:29:30  arno
% multiple condition baseline now in log space
%
% Revision 1.29  2007/05/19 01:15:09  toby
% error message incorrectly claiming precomputed data not found corrected
%
% Revision 1.28  2007/05/18 16:23:29  peter
% typo
%
% Revision 1.27  2007/05/18 16:17:54  peter
% typo
%
% Revision 1.26  2007/05/18 16:15:46  peter
% typo
%
% Revision 1.25  2007/05/18 16:07:09  peter
% added test for reading ERSP component vars.
%
% Revision 1.24  2007/01/26 18:06:35  arno
% nothing
%
% Revision 1.23  2006/11/10 02:28:46  arno
% bootstrap concatenation
%
% Revision 1.22  2006/11/03 19:45:06  arno
% do not mask ERSP by default
%
% Revision 1.21  2006/09/12 18:56:48  arno
% channel compatibility and many more features
%
% Revision 1.20  2006/05/13 17:46:22  arno
% transposing baseline to prevent crash
%
% Revision 1.19  2006/05/03 18:20:40  arno
% allowing to read data channels
%
% Revision 1.18  2006/03/29 17:47:10  scott
% help msg
%
% Revision 1.17  2006/03/28 14:58:10  arno
% reading ersp channel
%
% Revision 1.16  2006/03/22 00:46:32  scott
% help msg format only
%
% Revision 1.15  2006/03/14 03:02:51  scott
% help msg
%
% Revision 1.14  2006/03/14 02:23:18  scott
% help msg
%
% Revision 1.13  2006/03/14 01:59:38  scott
% help msg
%
% Revision 1.12  2006/03/13 19:09:03  arno
% no time range and freq range
%
% Revision 1.11  2006/03/11 07:28:49  arno
% header info
%
% Revision 1.10  2006/03/11 07:21:08  arno
% header
%
% Revision 1.9  2006/03/10 17:44:25  arno
% typo
%
% Revision 1.8  2006/03/10 15:49:19  arno
% fix reading ERSP
%
% Revision 1.7  2006/03/10 03:25:32  scott
% help msg -- ARNO, please check  -sm
%
% Revision 1.6  2006/03/10 00:39:39  arno
% error msg
%
% Revision 1.5  2006/03/09 23:29:34  arno
% implement new ERSP from Matlab and different structure ec...
%
% Revision 1.4  2006/03/09 19:37:06  arno
% header
%
% Revision 1.3  2006/03/08 20:34:24  arno
% rename func
%
% Revision 1.2  2006/03/07 22:16:23  arno
% use fullfile
%

function [STUDY, erspdata, alltimes, allfreqs] = std_readersp(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readersp;
    return;
end
if ~isstruct(ALLEEG) % old calling format
    [STUDY, erspdata] = std_readerspsub(STUDY, ALLEEG, varargin{:});
    return;
end;

STUDY = pop_erspparams(STUDY, 'default');
[opt moreopts] = finputcheck( varargin, { ...
    'type'          { 'string' 'cell' } { [] [] } '';
    'channels'      'cell'    []             {};
    'clusters'      'integer' []             [];
    'freqrange'     'real'    []             STUDY.etc.erspparams.freqrange;
    'timerange'     'real'    []             STUDY.etc.erspparams.timerange;
    'rmsubjmean'    'string'  { 'on' 'off' } 'off';
    'singletrials'  'string'  { 'on' 'off' } 'off';
    'subbaseline'   'string'  { 'on' 'off' }  STUDY.etc.erspparams.subbaseline;
    'component'     'integer' []               [];
    'infotype'      'string'  { 'ersp' 'itc' } 'ersp'; ...
    'datatype'      'string'  { 'ersp' 'itc' } 'ersp'; ...
    'subject'       'string'  []               '' }, ...
    'std_readersp', 'ignore');
if isstr(opt), error(opt); end;
nc = max(length(STUDY.condition),1);
ng = max(length(STUDY.group),1);
if ~strcmpi(opt.infotype, 'ersp'), opt.datatype = opt.infotype; end;
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
    eqtf  = isequal( STUDY.etc.erspparams.timerange, opt.timerange) && isequal( STUDY.etc.erspparams.freqrange, opt.freqrange);
    eqtfb = eqtf && isequal( STUDY.etc.erspparams.subbaseline, opt.subbaseline);
    if strcmpi(opt.singletrials,'off')
        if isfield(tmpstruct, [ dtype 'data']) && eqtfb && ~isempty(getfield(tmpstruct, [ dtype 'data']))
            dataread = 1;
        end;
    else
        if isfield(tmpstruct, [ dtype 'datatrials']) && eqtf && ~isempty(getfield(tmpstruct, [ dtype 'datatrials']))
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
        ersp     = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) );
        erspbase = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) );
        erspinds = cell( max(length(STUDY.condition),1), max(length(STUDY.group),1) );
        
        % find total nb of trials
        % -----------------------
        if strcmpi(opt.singletrials, 'on')
            tottrials = cell( length(STUDY.condition), length(STUDY.group) );
            for index = 1:length( STUDY.datasetinfo )
                condind = strmatch( STUDY.datasetinfo(index).condition, STUDY.condition );
                grpind  = strmatch( STUDY.datasetinfo(index).group    , STUDY.group     );
                if isempty(tottrials{condind, grpind}), tottrials{condind, grpind} = ALLEEG(index).trials;
                else       tottrials{condind, grpind} = tottrials{condind, grpind} + ALLEEG(index).trials;
                end;
            end;
        end;

        % read the data and select channels
        % ---------------------------------
        fprintf('Reading all %s data:', upper(dtype));
        for c = 1:nc
            for g = 1:ng
                %try     % this 'try' is a poor solution the problem of attempting
                    % to read specific channel/component data that doesn't exist
                    % called below by: allinds{c,g}(indtmp)
                    if strcmpi(opt.singletrials, 'on')
                        if ~isempty(opt.channels)
                            allsubjects = { STUDY.datasetinfo(:).subject };
                            if ~isempty(opt.subject) inds = strmatch( opt.subject, allsubjects(setinds{c,g}));
                            else inds = 1:length(allinds{c,g}); end;
                        else
                            if ~isempty(opt.component) inds = find( allinds{c,g} == STUDY.cluster(finalinds(ind)).comps(opt.component));
                            else inds = 1:length(allinds{c,g}); end;
                        end;
                        if ~isempty(inds)
                            count{c, g} = 1;
                            for indtmp = 1:length(inds)
                                [ tmpersp allfreqs alltimes tmpparams] = std_readtimef( ALLEEG, setinds{c,g}(inds(indtmp)), allinds{c,g}(inds(indtmp)), ...
                                    opt.timerange, opt.freqrange);
                                indices = [count{c, g}:count{c, g}+size(tmpersp,3)-1];
                                if indtmp == 1
                                     ersp{c, g} = permute(tmpersp, [2 1 3]);
                                else ersp{c, g}(:,:,indices) = permute(tmpersp, [2 1 3]);
                                end;
                                erspinds{c, g}(1:2,indtmp) = [ count{c, g} count{c, g}+size(tmpersp,3)-1 ];
                                count{c, g} = count{c, g}+size(tmpersp,3);
                                if size(tmpersp,3) ~= ALLEEG(setinds{c,g}(inds(indtmp))).trials
                                    error( sprintf('Wrong number of trials in datafile for dataset %d\n', setinds{c,g}(inds(indtmp))));
                                end;
                            end;
                        end;
                    elseif strcmpi(dtype, 'itc')
                        for indtmp = 1:length(allinds{c,g})
                            [ tmpersp allfreqs alltimes tmpparams] = std_readitc( ALLEEG, setinds{c,g}(indtmp), allinds{c,g}(indtmp), ...
                                opt.timerange, opt.freqrange);
                            if indtmp == 1, ersp{c, g} = repmat(single(0), [length(alltimes), length(allfreqs), length(allinds{c,g}) ]); end;
                            ersp{c, g}(:,:,indtmp)     = abs(permute(tmpersp    , [2 1]));
                        end;
                    else
                        for indtmp = 1:length(allinds{c,g})
                            [ tmpersp allfreqs alltimes tmpparams tmperspbase] = std_readerspsub( ALLEEG, setinds{c,g}(indtmp), allinds{c,g}(indtmp), ...
                                opt.timerange, opt.freqrange);
                            if indtmp == 1, ersp{c, g} = repmat(single(0), [length(alltimes), length(allfreqs), length(allinds{c,g}) ]); 
                                        erspbase{c, g} = repmat(single(0), [               1, length(allfreqs), length(allinds{c,g}) ]); 
                            end;
                            ersp{c, g}(    :,:,indtmp) = permute(tmpersp    , [2 1]);
                            erspbase{c, g}(:,:,indtmp) = 10*log(permute(tmperspbase, [2 1]));
                        end;
                    end
                %catch
                %end
                fprintf('.');
            end;
        end;
        fprintf('\n');

        % compute ERSP or ITC if trial mode
        % (since only the timef have been loaded)
        % ---------------------------------------
        if strcmpi(opt.singletrials, 'on')
            precomp.times = alltimes;
            precomp.freqs = allfreqs;
            precomp.recompute = dtype;
            for c = 1:nc
                for g = 1:ng
                    precomp.tfdata = permute(ersp{c,g}, [2 1 3]);                    
                    if strcmpi(dtype, 'itc')
                        [tmp ersp{c,g}] = newtimef(zeros(ALLEEG(1).pnts,2), ALLEEG(1).pnts, [ALLEEG(1).xmin ALLEEG(1).xmax]*1000, ...
                            ALLEEG(1).srate, tmpparams{:}, 'precomputed', precomp, 'verbose', 'off');
                    elseif strcmpi(dtype, 'ersp')
                        ersp{c,g} = newtimef(zeros(ALLEEG(1).pnts,2), ALLEEG(1).pnts, [ALLEEG(1).xmin ALLEEG(1).xmax]*1000, ...
                            ALLEEG(1).srate, tmpparams{:}, 'precomputed', precomp, 'verbose', 'off');
                    end;
                    ersp{c,g} = permute(ersp{c,g}, [2 1 3]);
                end;
            end;
        end;

        % compute average baseline across groups and conditions
        % -----------------------------------------------------
        if strcmpi(opt.subbaseline, 'on') && strcmpi(dtype, 'ersp') && strcmpi(opt.singletrials, 'off')
            disp('Recomputing baseline...');
            for g = 1:ng        % ng = number of groups
                for c = 1:nc    % nc = number of components
                    if strcmpi(opt.singletrials, 'on')
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
                    if strcmpi(opt.singletrials, 'on'), tmpmeanpowbase = repmat(meanpowbase, [length(alltimes) 1 tottrials{c,g}]);
                    else                                tmpmeanpowbase = repmat(meanpowbase, [length(alltimes) 1 1]);
                    end;
                    ersp{c,g} = ersp{c,g} - repmat(abs(erspbase{c,g}), [length(alltimes) 1 1 1]) + tmpmeanpowbase;
                end;
                clear meanpowbase;
            end;
        end;

%         if strcmpi(opt.statmode, 'common')
%             % collapse the two last dimensions before computing significance
%             % i.e. 18 subject with 4 channels -> same as 4*18 subjects
%             % --------------------------------------------------------------
%             disp('Using all channels for statistics...');
%             for c = 1:nc
%                 for g = 1:ng
%                     ersp{c,g} = reshape( ersp{c,g}, size(ersp{c,g},1), size(ersp{c,g},2), size(ersp{c,g},3)*size(ersp{c,g},4));
%                 end;
%             end;
%         end;

        % copy data to structure
        % ----------------------
        tmpstruct = setfield(tmpstruct, [ dtype 'freqs' ], allfreqs);
        tmpstruct = setfield(tmpstruct, [ dtype 'times' ], alltimes);
        if strcmpi(opt.singletrials, 'on')
            tmpstruct = setfield(tmpstruct, [ dtype 'datatrials' ], ersp);
            tmpstruct = setfield(tmpstruct, [ dtype 'subjinds' ], erspinds);
            tmpstruct = setfield(tmpstruct, [ dtype 'times' ], alltimes);
            if ~isempty(opt.channels)
                 tmpstruct = setfield(tmpstruct, [ dtype 'trialinfo' ], opt.subject);
            else tmpstruct = setfield(tmpstruct, [ dtype 'trialinfo' ], opt.component);
            end;
        else
            tmpstruct = setfield(tmpstruct, [ dtype 'data' ], ersp);
            if strcmpi(dtype, 'ersp')
                tmpstruct = setfield(tmpstruct, [ dtype 'base' ], erspbase);
            end;
        end;
        
        % copy results to structure
        % -------------------------
        fieldnames = { [ dtype 'data' ] [ dtype 'freqs' ] [ dtype 'datatrials' ] [ dtype 'subjinds' ] ...
                       [ dtype 'base' ] [ dtype 'times' ] [ dtype 'trialinfo' ] 'allinds' 'setinds' };
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
    erspdata = cell(nc, ng);
    for ind =  1:length(erspdata(:))
        if strcmpi(opt.singletrials, 'on')
             tmpdat = getfield(structdat(allinds(1)), [ dtype 'datatrials' ]);
        else tmpdat = getfield(structdat(allinds(1)), [ dtype 'data' ]);
        end;
        if ndims(tmpdat{ind}) == 2, erspdata{ind} = zeros([ size(tmpdat{ind}) 1 length(allinds)]);
        else                        erspdata{ind} = zeros([ size(tmpdat{ind}) length(allinds)]);
        end;
        for index = 1:length(allinds)
            if strcmpi(opt.singletrials, 'on')
                 tmpdat = getfield(structdat(allinds(index)), [ dtype 'datatrials' ]);
            else tmpdat = getfield(structdat(allinds(index)), [ dtype 'data' ]);
            end;
            erspdata{ind}(:,:,:,index) = tmpdat{ind};
            allfreqs = getfield(structdat(allinds(index)), [ dtype 'freqs' ]);
            alltimes = getfield(structdat(allinds(index)), [ dtype 'times' ]);
            compinds = structdat(allinds(index)).allinds;
            setinds  = structdat(allinds(index)).setinds;
        end;
        erspdata{ind} = permute(erspdata{ind}, [1 2 4 3]); % time freqs elec subjects
    end;
    if ~isempty(opt.subject) && strcmpi(opt.singletrials,'off')
        erspdata = std_selsubject(erspdata, opt.subject, setinds, { STUDY.datasetinfo(:).subject }, 2); 
    end;
else
    if strcmpi(opt.singletrials, 'on')
         erspdata = getfield(STUDY.cluster(allinds(1)), [ dtype 'datatrials' ]);
    else erspdata = getfield(STUDY.cluster(allinds(1)), [ dtype 'data' ]);
    end;
    allfreqs = getfield(STUDY.cluster(allinds(1)), [ dtype 'freqs' ]);
    alltimes = getfield(STUDY.cluster(allinds(1)), [ dtype 'times' ]);
    compinds = STUDY.cluster(allinds(1)).allinds;
    setinds  = STUDY.cluster(allinds(1)).setinds;
    if ~isempty(opt.component) && length(allinds) == 1 && strcmpi(opt.singletrials,'off')
        erspdata = std_selcomp(STUDY, erspdata, allinds, setinds, compinds, opt.component);
    end;
end;

% std_readerspsub() - Returns the equal-log-frequency spaced mean event-related spectral 
%                  perturbation(s) (ERSP(s)) for a requested ICA component. The ERSP 
%                  is assumed to have been saved in a Matlab file, 
%                          [dataset_name].icaersp 
%                  in the same folder as the dataset file.
%                  If this file does not exist, use std_ersp() to create it, else 
%                  use a pre-clustering function: pop_preclust() or std_preclust(), 
%                  that calls it. Interpretation of the ERSP requires some input 
%                  variables used to compute it: frequency range, window width, 
%                  resolution, probability threshold, and wavelet type (FFT or 
%                  wavelet_cycles). See >> timef help and >> timef details
% Usage:    
%     >>  [logersp, logfreqs, timevals, params, baseersp] = std_readersp(ALLEEG, setindex, component, ...
%                                                            time_range, freq_range);  
% Inputs:
%   ALLEEG     - vector of EEG datasets (can also be one EEG dataset). Must contain 
%                the dataset of interest (see 'setind' below).
%   setindex   - [integer] index of the EEG dataset in the ALLEEG structure for 
%                which to read a component log-frequency ERSP.
%   component  - [integer] component index in the selected EEG dataset for which 
%                to read the ERSP 
%   time_range - [min max ms] ERSP time (latency) range of interest
%   freq_range - [min max Hz] ERSP frequency range of interest
%
% Outputs:
%   logersp    - the log-frequency ERSP for the requested ICA component 
%                in the specified dataset. Dimensions: (equal log-spaced) 
%                frequencies by epoch latencies (unit: dB diff from baseline)
%   logfreqs   - vector of equal-log-spaced ERSP frequencies (Hz) 
%   timevals      - vector of ERSP times (latencies) (s)
%   params     - structure of timef() parameters saved with the ERSP
%   baseersp   - condition-specific baseline. Third dimesnsion corresponds to 
%                conditions.
%   
% See also:  std_ersp(), pop_preclust(), std_preclust(), timef()
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, February, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004, arno@sccn.ucsd.edu
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
% Revision 1.36  2010/02/16 08:43:21  arno
% New single-trial reading/writing
%

function [logersp, logfreqs, timevals, params, baseersp] = std_readerspsub(ALLEEG, abset, comp, timewindow, freqrange);

if nargin < 4
    timewindow = [];
end;
if nargin < 5
    freqrange = [];
end;

% multiple entry
% --------------
if (length(comp) > 1 & comp(1) > 0) | length(comp) > length(abset) % recursive call if multiple components
    for index = 1:length(comp)
        [tmpersp, logfreqs, timevals, params, tmpbase] = std_readerspsub(ALLEEG, abset, comp(index), timewindow, freqrange);
        logersp(index,:,:,:) = tmpersp;
        baseersp(index,:,:)  = tmpbase;
    end;
    return;
end;

for k = 1: length(abset)    
    
    if comp < 0
        filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'datersp']);
        comptmp  = -comp(k);
        prefix   = 'chan';
    else    
        filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'icaersp']);
        prefix   = 'comp';
        comptmp  = comp;
    end;
    try
        tmpersp   = load( '-mat', filename, 'parameters', 'times', 'freqs');
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end;
    
    tmpersp.parameters = removedup(tmpersp.parameters);
    params    = struct(tmpersp.parameters{:});
    params.times = tmpersp.times;
    params.freqs = tmpersp.freqs;
    tlen         = length(tmpersp.times);
    flen         = length(tmpersp.freqs);
    if isempty(comp)
        logersp   = [];
        logfreqs  = [];
        timevals  = [];
        return;
    end;
    if isinf(comp) % only read time and freqs
        abset = [];
        erspallboot{1} = [];
        logersp = [];
    else
        tmpersp2   = load( '-mat', filename, ...
                         [ prefix int2str(comptmp) '_ersp'], ...
                         [ prefix int2str(comptmp) '_erspbase'], ...
                         [ prefix int2str(comptmp) '_erspboot']);
        if ~isfield(tmpersp2,[ prefix int2str(comptmp) '_ersp']) |   ...
           ~isfield(tmpersp2,[ prefix int2str(comptmp) '_erspbase']) |   ...
           ~isfield(tmpersp2,[ prefix int2str(comptmp) '_erspboot']) 
          fprintf('\nERSP data for component %d not found!\n       In file %s\n',comptmp,filename); 
        end
        erspall{k}     = double(getfield(tmpersp2, [ prefix int2str(comptmp) '_ersp']));
        erspallboot{k} = double(getfield(tmpersp2, [ prefix int2str(comptmp) '_erspboot']));
        erspallbase{k} = double(getfield(tmpersp2, [ prefix int2str(comptmp) '_erspbase']));
    end;
end

% compute average baseline across conditions
% ------------------------------------------
if length(abset) > 1 
    % mean baseline for requested component across conditions
    % -------------------------------------------------------
	ave_baseline = zeros(size(erspallbase{1}(:)')); 
	for cond = 1:length(abset)
        ave_baseline = ave_baseline + erspallbase{cond}(:)'/length(abset);
	end    
    
    % apply mean baseline
    % -------------------
    for cond = 1:length(abset)
        % add back former baseline to the ERSP and subtract new baseline
        if any(erspallbase{cond}) > 100, 
            error([ 'You must recompute time-frequency decomposition' 10 ...
                    'changes greater than 100dB detected in ERSP which means' 10 ...
                    'that you have most likely computed these ERSP using a' 10 ...
                    'version of EEGLAB < 6.00' ]);
        end;
        
        erspall{cond} = erspall{cond} + repmat( erspallbase{cond}(:),[1 tlen]);
        erspall{cond} = erspall{cond} - repmat( ave_baseline(:)     ,[1 tlen]);
        
        % same for bootstrap array
        if ~isempty(erspallboot{cond})
            erspallboot{cond} = erspallboot{cond} + repmat(erspallbase{cond}',[1 2]);
            erspallboot{cond} = erspallboot{cond} - repmat(ave_baseline'     ,[1 2]);  
        end;
    end
end

% select plotting or clustering time/freq range
% ---------------------------------------------
if ~isempty(timewindow)
    if timewindow(1) > tmpersp.times(1) | timewindow(end) < tmpersp.times(end)
        maxind = max(find(tmpersp.times <= timewindow(end)));
        minind = min(find(tmpersp.times >= timewindow(1)));
    else
        minind = 1;
        maxind = tlen;
    end
else
    minind = 1;
    maxind = tlen;
end
if ~isempty(freqrange)
    if freqrange(1) > exp(1)^tmpersp.freqs(1) | freqrange(end) < exp(1)^tmpersp.freqs(end)
        fmaxind = max(find(tmpersp.freqs <= freqrange(end)));
        fminind = min(find(tmpersp.freqs >= freqrange(1)));
    else
        fminind = 1;
        fmaxind = flen;
    end
else
    fminind = 1;
    fmaxind = flen;
end

% Mask ERSP
% ---------
%if ~isempty(erspallboot{1})
%    for cond  = 1:length(abset)
%        minersp= repmat(erspallboot{cond}(1,:)',1,size(erspall{1},2));
%        maxersp= repmat(erspallboot{cond}(2,:)',1,size(erspall{1},2));
%        erspall{cond}(find(erspall{cond}<maxersp & erspall{cond}>minersp)) = 0;
%    end
%end;

% return parameters
% ----------------
for cond  = 1:length(abset)
    ersp = erspall{cond}(fminind:fmaxind,minind:maxind);
    logersp(:,:,cond) = ersp;
    baseersp(:,cond)  = erspallbase{cond}(fminind:fmaxind)';
end;
logfreqs = tmpersp.freqs(fminind:fmaxind);
timevals = tmpersp.times(minind:maxind);

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
    [tmp indices] = unique(cella(1:2:end));
    if length(tmp) ~= length(cella)/2
        %fprintf('Warning: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
    end;
    cella = cella(sort(union(indices*2-1, indices*2)));
    
% std_readtimef() - returns the log-frequency time-frequency decomposition for a 
%                 specified ICA component. The component time-freqs for the dataset 
%                 are assumed to have been saved in a Matlab file, 
%                 [dataset_)name].icatimef, in the same directory as the dataset.
%                 If no such file exists, use std_ersp() to create it, else 
%                 the pre-clustering functions that call it: pop_preclust, 
%                 std_preclust().  The input variables used to compute the
%                 time-freq. are returned: frequency_range, time_range, resolution, 
%                 probability_threshold, and wavelet type (FFT | wavelet 
%                 cycles). See timef() for details. 
% Usage:    
%   >> [logtimef, logfreqs, times] = std_readtimef(ALLEEG, setindx, component, ...
%                                                       time_range, freq_range);  
% Inputs:
%   ALLEEG     - EEG dataset vector (can also be an EEG set). 
%                Must contain the dataset of interest (see 'setindx' below).
%   setindx    -  [integer] index of the EEG dataset in ALLEEG for which 
%                 to return the log-frequency time-freq decomposition.
%   component  - [integer] component index in the selected EEG dataset for 
%                which to return the time-frequency decomposition. 
%   time_range - [min max in ms] time window 
%   freq_range - [min max in Hz] frequency range 
%
% Outputs:
%   logtimef   - the equal log-spaced frequency time-frequency decomposition
%                for the requested ICA component in the specified dataset. 
%                Its dimensions are (equal log-spaced) frequencies by times. 
%   logfreqs   - vector of log-spaced frequencies, in Hz
%   times      - vector of times (latencies), in ms.
%   params     - full structure of time-freq. parameters saved
%
%  See also  std_ersp(), std_readersp(), pop_preclust(), eeg_preclust(), 
%               eeg_createdata()
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, May, 2006

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
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
% Revision 1.36  2010/02/16 08:43:21  arno
% New single-trial reading/writing
%

function [logtimef, logfreqs, timevals, params] = std_readtimef(ALLEEG, abset, comp, timewindow, freqrange);

if nargin < 4
    timewindow = [];
end;
if nargin < 5
    freqrange = [];
end;

% multiple entry
% --------------
if length(comp) > 1
    for index = 1:length(comp)
        [tmptimef, logfreqs, timevals, params] = std_readtimef(ALLEEG, abset, comp(index), timewindow, freqrange);
        logtimef(index,:,:,:) = tmptimef;
    end;
    return;
end;

for k = 1: length(abset)    
    
    if comp < 0
        filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'dattimef']);
        comp   = -comp;
        prefix = 'chan';
    else    
        filename = fullfile( ALLEEG(abset(k)).filepath,[ ALLEEG(abset(k)).filename(1:end-3) 'icatimef']);
        prefix = 'comp';
    end;
    try
        tmpersp   = load( '-mat', filename, 'parameters', 'times', 'freqs');
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end;
    params    = tmpersp.parameters;
    if isempty(comp)
        logtimef    = [];
        logfreqs  = [];
        timevals  = [];
        return;
    end;
    tmptimef  = load( '-mat', filename, 'parameters', 'times', 'freqs', ...
                     [ prefix int2str(comp) '_timef']);
    
    tlen        = length(tmptimef.times);
    flen        = length(tmptimef.freqs);
    timefall{k} = double(getfield(tmptimef, [ prefix int2str(comp) '_timef']));
    
end

% select plotting or clustering time/freq range
% ---------------------------------------------
if ~isempty(timewindow)
    if timewindow(1) > tmptimef.times(1) | timewindow(end) < tmptimef.times(end)
        maxind = max(find(tmptimef.times <= timewindow(end)));
        minind = min(find(tmptimef.times >= timewindow(1)));
    else
        minind = 1;
        maxind = tlen;
    end
else
    minind = 1;
    maxind = tlen;
end
if ~isempty(freqrange)
    if freqrange(1) > exp(1)^tmptimef.freqs(1) | freqrange(end) < exp(1)^tmptimef.freqs(end)
        fmaxind = max(find(tmptimef.freqs <= freqrange(end)));
        fminind = min(find(tmptimef.freqs >= freqrange(1)));
    else
        fminind = 1;
        fmaxind = flen;
    end
else
    fminind = 1;
    fmaxind = flen;
end

% return parameters
% ----------------
for cond  = 1:length(abset)
    logtimef(:,:,:,cond) = timefall{cond}(fminind:fmaxind,minind:maxind,:);
end;
logfreqs = tmptimef.freqs(fminind:fmaxind);
timevals = tmptimef.times(minind:maxind);
