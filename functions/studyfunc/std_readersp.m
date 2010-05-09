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

% $Log: std_readersp.m,v $
% Revision 1.38  2010/03/09 06:19:45  arno
% fixed reading single trials for multiple subjects
%
% Revision 1.37  2010/02/24 10:52:37  arno
% Implemented new single trial statistics
%
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
    'design'        'integer' []             STUDY.currentdesign;    
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
if ~strcmpi(opt.infotype, 'ersp'), opt.datatype = opt.infotype; end;
nc = max(length(STUDY.design(opt.design).condition),1);
ng = max(length(STUDY.design(opt.design).group),1);
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
        setinds       = tmpstruct.setinds;
        for i=1:length(allinds(:)), allinds{i} = -allinds{i}; end; % invert sign for reading
    else
        tmpstruct = STUDY.cluster(finalinds(ind));
        allinds       = tmpstruct.allinds;
        setinds       = tmpstruct.setinds;
    end;

    % check if data is already here
    % -----------------------------
    dataread = 0;
    eqtf  = isequal( STUDY.etc.erspparams.timerange, opt.timerange) && ...
            isequal( STUDY.etc.erspparams.freqrange, opt.freqrange);
    eqtfb = isequal( STUDY.etc.erspparams.subbaseline, opt.subbaseline) && eqtf;
    if strcmpi(opt.singletrials,'off')
        if isfield(tmpstruct, [ dtype 'data']) && eqtfb && ~isempty(getfield(tmpstruct, [ dtype 'data']))
            dataread = 1;
        end;
    else
        if isfield(tmpstruct, [ dtype 'datatrials']) && eqtf && ~isempty(getfield(tmpstruct, [ dtype 'datatrials']))
            if ~isempty(opt.channels) && strcmpi(getfield(tmpstruct, [ dtype 'trialinfo' ]), opt.subject)
                dataread = 1; 
            elseif isempty(opt.channels) && isequal(getfield(tmpstruct, [ dtype 'trialinfo' ]), opt.component) 
                dataread = 1; 
            end;
        end;
    end;
    
    if ~dataread
        
        % reserve arrays
        % --------------
        ersp     = cell( nc, ng );
        erspbase = cell( nc, ng );
        erspinds = cell( nc, ng );
        
        % find total nb of trials
        % -----------------------
        setinfo = STUDY.design(opt.design).setinfo;
        if strcmpi(opt.singletrials, 'on')
            tottrials = cell( nc, ng );
            for index = 1:length(STUDY.design(opt.design).setinfo)
                condind = strmatch( setinfo(index).condition, STUDY.design(opt.design).condition );
                grpind  = strmatch( setinfo(index).group    , STUDY.design(opt.design).group     );
                if isempty(tottrials{condind, grpind}), tottrials{condind, grpind} = ALLEEG(index).trials;
                else       tottrials{condind, grpind} = tottrials{condind, grpind} + ALLEEG(index).trials;
                end;
            end;
        end;

        % read the data and select channels
        % ---------------------------------
        fprintf('Reading all %s data:', upper(dtype));
        setinfo = STUDY.design(opt.design).setinfo;        
        %try     % this 'try' is a poor solution the problem of attempting
        % to read specific channel/component data that doesn't exist
        % called below by: allinds{c,g}(indtmp)
        if strcmpi(opt.singletrials, 'on')
            for c = 1:nc
                for g = 1:ng
                    if ~isempty(opt.channels)
                        allsubjects = { STUDY.design(opt.design).setinfo.subject };
                        if ~isempty(opt.subject), inds = strmatch( opt.subject, allsubjects(setinds{c,g}));
                        else inds = 1:length(allinds{c,g}); end;
                    else
                        if ~isempty(opt.component) inds = find( allinds{c,g} == STUDY.cluster(finalinds(ind)).comps(opt.component));
                        else inds = 1:length(allinds{c,g}); end;
                    end;
                    if ~isempty(inds)
                        count{c, g} = 1;
                        for indtmp = 1:length(inds)
                            setindtmp = STUDY.design(opt.design).setinfo(setinds{c,g}(inds(indtmp))).setindex;
                            [ tmpersp tmpparams alltimes allfreqs] = std_readfile(setindtmp, 'measure', 'timef', 'dataindices', allinds{c,g}(inds(indtmp)), 'timelimits', opt.timerange, 'freqlimits', opt.freqrange);
                            indices = [count{c, g}:count{c, g}+size(tmpersp,3)-1];
                            if indtmp == 1
                                 ersp{c, g} = permute(tmpersp, [2 1 3]);
                            else ersp{c, g}(:,:,indices) = permute(tmpersp, [2 1 3]);
                            end;
                            erspinds{c, g}(1:2,indtmp) = [ count{c, g} count{c, g}+size(tmpersp,3)-1 ];
                            count{c, g} = count{c, g}+size(tmpersp,3);
                            if size(tmpersp,3) ~= ALLEEG(setindtmp).trials
                                %error( sprintf('Wrong number of trials in datafile for dataset %d\n', setinds{c,g}(inds(indtmp))));
                            end;
                        end;
                    end;
                end;
            end;
        else
            for c = 1:nc
                for g = 1:ng
                    options = { 'dataindices', allinds{c,g}(:), 'timelimits', opt.timerange, 'freqlimits', opt.freqrange };
                    if strcmpi(dtype, 'ersp') 
                         erspbase{c, g}                             = std_readfile( setinfo(setinds{c,g}(:)), 'measure', 'erspbase', options{:});
                         [ ersp{c, g} tmpparams alltimes allfreqs ] = std_readfile( setinfo(setinds{c,g}(:)), 'measure', 'ersp'    , options{:});
                    else [ ersp{c, g} tmpparams alltimes allfreqs ] = std_readfile( setinfo(setinds{c,g}(:)), 'measure', 'itc'     , options{:});
                         ersp{c, g} = abs(ersp{c, g});
                    end;
                    fprintf('.');
                    %ersp{c, g}      = permute(ersp{c, g}           , [3 2 1]);
                    %erspbase{c, g}  = 10*log(permute(erspbase{c, g}, [3 2 1]));
                end;
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
                    if ~isempty(ersp{c,g})
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
                erspbasetmp = reshape(erspbase{c,g}, [size(meanpowbase,1) 1 size(meanpowbase,2)]);
                meanpowbase = reshape(meanpowbase  , [size(meanpowbase,1) 1 size(meanpowbase,2)]);
                
                % subtract average baseline
                % -------------------------
                for c = 1:nc
                    if strcmpi(opt.singletrials, 'on'), tmpmeanpowbase = repmat(meanpowbase, [1 length(alltimes) tottrials{c,g}]);
                    else                                tmpmeanpowbase = repmat(meanpowbase, [1 length(alltimes) 1]);
                    end;
                    ersp{c,g} = ersp{c,g} - repmat(abs(erspbasetmp), [1 length(alltimes) 1 1]) + tmpmeanpowbase;
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
