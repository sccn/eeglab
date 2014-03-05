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
%  'subject'   - [string] select a specific subject {default:all}
%  'component' - [integer] select a specific component in a cluster
%                 {default:all}
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

function [STUDY, erspdata, alltimes, allfreqs, erspbase, events, unitPower] = std_readersp(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_readersp;
    return;
end
events = {};
unitPower = 'dB';
if ~isstruct(ALLEEG) % old calling format
    dataset = ALLEEG;
    EEG = STUDY(dataset);
    comp = varargin{1};
    if length(varargin) > 1, timelim = varargin{2}; else timelim = []; end;
    if length(varargin) > 2, freqlim = varargin{3}; else freqlim = []; end;
    filebase = fullfile(EEG.filepath, EEG.filename(1:end-4));
    if comp < 0
        error('Old format function call, channel reading not supported');
        [ersp params alltimes allfreqs] = std_readfile(filebase, 'channels', -comp, 'timelimits', timelim, 'freqlimits', freqlim);
    else
        [ersp     params alltimes allfreqs] = std_readfile(filebase, 'components', comp, 'timelimits', timelim, 'freqlimits', freqlim, 'measure', 'ersp');
        [erspbase params alltimes allfreqs] = std_readfile(filebase, 'components', comp, 'timelimits', timelim, 'freqlimits', freqlim, 'measure', 'erspbase');
    end;
    STUDY = ersp;
    erspdata = allfreqs;
    allfreqs = params;
    return;
end;

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

for ind = 1:length(finalinds)
    erspbase = cell( nc, ng );

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
    if strcmpi(dtype, 'erpim') % ERP images
        dataread = 0;
        eqtf  = isequal( STUDY.etc.erpimparams.timerange , opt.timerange) && ...
                isequal( STUDY.etc.erpimparams.trialrange, opt.trialrange);
        if isfield(tmpstruct, 'erpimdata') && eqtf && ~isempty(tmpstruct.erpimdata)
            dataread = 1;
        end;
    else
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
    end;
    
    if ~dataread
        
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
            ersp     = cell( nc, ng );
            erspinds = cell( nc, ng );
            
            % find total nb of trials
            % THIS CODE IS NOT NECESSARY ANY MORE (SEE BUG 1170)
            % -----------------------
            %setinfo = STUDY.design(opt.design).cell;
            %tottrials = cell( nc, ng );
            %if strcmpi(opt.singletrials, 'on')
            %    for indSet = 1:length(setinfo)
            %        condind = std_indvarmatch( setinfo(indSet).value{1}, STUDY.design(opt.design).variable(1).value );
            %        grpind  = std_indvarmatch( setinfo(indSet).value{2}, STUDY.design(opt.design).variable(2).value );
            %        if isempty(tottrials{condind, grpind}), tottrials{condind, grpind} = sum(cellfun(@length, setinfo(indSet).trials));
            %        else       tottrials{condind, grpind} = tottrials{condind, grpind} + sum(cellfun(@length, setinfo(indSet).trials));
            %        end;
            %    end;
            %end;
            
            % read the data and select channels
            % ---------------------------------
            fprintf('Reading all %s data:', upper(dtype));
            setinfo = STUDY.design(opt.design).cell;
            if strcmpi(opt.singletrials, 'on')
                for c = 1:nc
                    for g = 1:ng
                        if ~isempty(opt.channels)
                            allsubjects = { STUDY.design(opt.design).cell.case };
                            if ~isempty(opt.subject), inds = strmatch( opt.subject, allsubjects(setinds{c,g}));
                            else inds = 1:length(allinds{c,g}); end;
                        else
                            if ~isempty(opt.component) inds = find( allinds{c,g} == STUDY.cluster(finalinds(ind)).comps(opt.component));
                            else inds = 1:length(allinds{c,g}); end;
                        end;
                        if ~isempty(inds)
                            count{c, g} = 1;
                            for indtmp = 1:length(inds)
                                setindtmp = STUDY.design(opt.design).cell(setinds{c,g}(inds(indtmp)));
                                tmpopts = { 'measure', 'timef' 'timelimits', opt.timerange, 'freqlimits', opt.freqrange };
                                if ~isempty(opt.channels), [ tmpersp tmpparams alltimes allfreqs] = std_readfile(setindtmp, 'channels', allChangrp(allinds{c,g}(inds(indtmp))), tmpopts{:});
                                else                       [ tmpersp tmpparams alltimes allfreqs] = std_readfile(setindtmp, 'components',          allinds{c,g}(inds(indtmp)),  tmpopts{:});
                                end;
                                indices = [count{c, g}:count{c, g}+size(tmpersp,3)-1];
                                if indtmp == 1
                                    ersp{c, g} = permute(tmpersp, [2 1 3]);
                                else ersp{c, g}(:,:,indices) = permute(tmpersp, [2 1 3]);
                                end;
                                erspinds{c, g}(1:2,indtmp) = [ count{c, g} count{c, g}+size(tmpersp,3)-1 ];
                                count{c, g} = count{c, g}+size(tmpersp,3);
                                if size(tmpersp,3) ~= sum(cellfun(@length, setindtmp.trials))
                                    error( sprintf('Wrong number of trials in datafile for design index %d\n', setinds{c,g}(inds(indtmp))));
                                end;
                            end;
                        end;
                    end;
                end;
            else
                for c = 1:nc
                    for g = 1:ng
                        if ~isempty(setinds{c,g})
                            if ~isempty(opt.channels), opts = { 'channels',  allChangrp(allinds{c,g}(:)), 'timelimits', opt.timerange, 'freqlimits', opt.freqrange };
                            else                       opts = { 'components',           allinds{c,g}(:) , 'timelimits', opt.timerange, 'freqlimits', opt.freqrange };
                            end;
                            if strcmpi(dtype, 'ersp')
                                 erspbase{c, g}                             = std_readfile( setinfo(setinds{c,g}(:)), 'measure', 'erspbase', opts{:});
                                 [ ersp{c, g} tmpparams alltimes allfreqs ] = std_readfile( setinfo(setinds{c,g}(:)), 'measure', 'ersp'    , opts{:});
                            else [ ersp{c, g} tmpparams alltimes allfreqs ] = std_readfile( setinfo(setinds{c,g}(:)), 'measure', 'itc'     , opts{:});
                                 ersp{c, g} = abs(ersp{c, g});
                            end;
                            fprintf('.');
                            %ersp{c, g}      = permute(ersp{c, g}           , [3 2 1]);
                            %erspbase{c, g}  = 10*log(permute(erspbase{c, g}, [3 2 1]));
                        end;
                    end;
                end;
            end;
            fprintf('\n');
            
            % compute ERSP or ITC if trial mode
            % (since only the timef have been loaded)
            % ---------------------------------------
            if strcmpi(opt.singletrials, 'on')
                tmpparams2      = fieldnames(tmpparams); tmpparams2 = tmpparams2';
                tmpparams2(2,:) = struct2cell(tmpparams);
                precomp.times = alltimes;
                precomp.freqs = allfreqs;
                precomp.recompute = dtype;
                for c = 1:nc
                    for g = 1:ng
                        if ~isempty(ersp{c,g})
                            precomp.tfdata = permute(ersp{c,g}, [2 1 3]);
                            if strcmpi(dtype, 'itc')
                                [tmp ersp{c,g}] = newtimef(zeros(ALLEEG(1).pnts,2), ALLEEG(1).pnts, [ALLEEG(1).xmin ALLEEG(1).xmax]*1000, ...
                                    ALLEEG(1).srate, [], tmpparams2{:}, 'precomputed', precomp, 'verbose', 'off');
                            elseif strcmpi(dtype, 'ersp')
                                [ersp{c,g} tmp] = newtimef(zeros(ALLEEG(1).pnts,2), ALLEEG(1).pnts, [ALLEEG(1).xmin ALLEEG(1).xmax]*1000, ...
                                    ALLEEG(1).srate, [], tmpparams2{:}, 'precomputed', precomp, 'verbose', 'off');
                            end;
                            ersp{c,g} = permute(ersp{c,g}, [2 1 3]);
                        end;
                    end;
                end;
            end;
            
            % compute average baseline across groups and conditions
            % -----------------------------------------------------
            if strcmpi(opt.subbaseline, 'on') && strcmpi(dtype, 'ersp')
                if strcmpi(opt.singletrials, 'on')
                    disp('WARNING: no ERSP baseline may not be subtracted when using single trials');
                else
                    disp('Recomputing baseline...');
                    if strcmpi(paired1, 'on') && strcmpi(paired2, 'on')
                        disp('Removing ERSP baseline for both indep. variables');
                        meanpowbase = computeerspbaseline(erspbase(:), opt.singletrials);
                        ersp        = removeerspbaseline(ersp, erspbase, meanpowbase);
                    elseif strcmpi(paired1, 'on')
                        disp('Removing ERSP baseline for first indep. variables (second indep. var. is unpaired)');
                        for g = 1:ng        % ng = number of groups
                            meanpowbase = computeerspbaseline(erspbase(:,g), opt.singletrials);
                            ersp(:,g)   = removeerspbaseline(ersp(:,g), erspbase(:,g), meanpowbase);
                        end;
                    elseif strcmpi(paired2, 'on')
                        disp('Removing ERSP baseline for second indep. variables (first indep. var. is unpaired)');
                        for c = 1:nc        % ng = number of groups
                            meanpowbase = computeerspbaseline(erspbase(c,:), opt.singletrials);
                            ersp(c,:)   = removeerspbaseline(ersp(c,:), erspbase(c,:), meanpowbase);
                        end;
                    else
                        disp('Not removing ERSP baseline (both indep. variables are unpaired');
                    end;
                end;
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
        if ~isempty(events)
            tmpstruct = setfield(tmpstruct, [ dtype 'events' ], events);
        end;
        tmpstruct = setfield(tmpstruct, [ dtype ordinate ], allfreqs);
        tmpstruct = setfield(tmpstruct, [ dtype 'times'  ], alltimes);
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
        fields = { [ dtype 'data'     ] [ dtype 'events' ] [ dtype ordinate ] [ dtype 'datatrials' ] ...
                   [ dtype 'subjinds' ] [ dtype 'base'   ] [ dtype 'times'  ] [ dtype 'trialinfo'  ]  'allinds' 'setinds' };
        for f = 1:length(fields)
            if isfield(tmpstruct, fields{f}),
                tmpdata = getfield(tmpstruct, fields{f});
                if ~isempty(opt.channels)
                     STUDY.changrp = setfield(STUDY.changrp, { finalinds(ind) }, fields{f}, tmpdata);
                else STUDY.cluster = setfield(STUDY.cluster, { finalinds(ind) }, fields{f}, tmpdata);
                end;
            end;
        end;
    end;
end;

% output unit
% -----------
if exist('tmpparams') ~= 1
    if ~isempty(opt.channels), [tmpersp tmpparams] = std_readfile(STUDY.design(opt.design).cell(1), 'channels', {STUDY.changrp(1).name}, 'measure', opt.datatype);
    else                       [tmpersp tmpparams] = std_readfile(STUDY.design(opt.design).cell(1), 'components', STUDY.cluster(finalinds(end)).allinds{1,1}(1), 'measure', opt.datatype);
    end;
end;
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
if ~isempty(opt.channels)
    structdat = STUDY.changrp;
    erspdata  = cell(nc, ng);
    events    = {};
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
            allfreqs = getfield(structdat(allinds(index)), [ dtype ordinate ]);
            alltimes = getfield(structdat(allinds(index)), [ dtype 'times'  ]);
            if isfield(structdat, [ dtype 'events' ])
                events   = getfield(structdat(allinds(index)), [ dtype 'events' ]);
            end;
            compinds = structdat(allinds(index)).allinds;
            setinds  = structdat(allinds(index)).setinds;
        end;
        erspdata{ind} = permute(erspdata{ind}, [1 2 4 3]); % time freqs elec subjects
    end;
    if ~isempty(opt.subject) && strcmpi(opt.singletrials,'off')
        erspdata = std_selsubject(erspdata, opt.subject, setinds, { STUDY.design(opt.design).cell.case }, 2); 
    end;
else
    if strcmpi(opt.singletrials, 'on')
         erspdata = getfield(STUDY.cluster(allinds(1)), [ dtype 'datatrials' ]);
    else erspdata = getfield(STUDY.cluster(allinds(1)), [ dtype 'data' ]);
    end;
    allfreqs = getfield(STUDY.cluster(allinds(1)), [ dtype ordinate ]);
    alltimes = getfield(STUDY.cluster(allinds(1)), [ dtype 'times'  ]);
    if isfield(STUDY.cluster, [ dtype 'events' ])
        events   = getfield(STUDY.cluster(allinds(1)), [ dtype 'events' ]);
    end;
    compinds = STUDY.cluster(allinds(1)).allinds;
    setinds  = STUDY.cluster(allinds(1)).setinds;
    if ~isempty(opt.component) && length(allinds) == 1 && strcmpi(opt.singletrials,'off')
        erspdata = std_selcomp(STUDY, erspdata, allinds, setinds, compinds, opt.component);
    end;
end;

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
   