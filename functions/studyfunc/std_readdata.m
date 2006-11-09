% std_readdata() - load one or more requested measures 
%                   ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map']
%                   for all components of a specified cluster.  
%                   Called by cluster plotting functions: std_envtopo(), 
%                   std_erpplot(), std_erspplot(), ...
% Usage:
%         >> [STUDY clustinfo] = std_clustread(STUDY, ALLEEG, cluster, infotype);
% Inputs:
%         STUDY - studyset structure containing some or all files in ALLEEG
%        ALLEEG - vector of loaded EEG datasets
%       cluster - cluster number in STUDY
%      infotype - ['erp'|'spec'|'ersp'|'itc'|'dipole'|'map'] type of stored
%                 cluster information to read. May also be a cell array of
%                 these types, for example: { 'erp' 'map' 'dipole' }.
%                 ************* 'map' and 'dipole' not implemented yet
%
% Output:
%    STUDY     - updated STUDY structure
%    clustinfo - structure of specified cluster information. This is 
%                equal to STUDY.cluster(cluster)     
%
%         clustinfo.erpdata   % (ncomps, ntimes) array of component ERPs
%         clustinfo.erptimes  % vector of ERP epoch latencies
%
%         clustinfo.specdata  % (ncomps, nfreqs) array of component spectra
%         clustinfo.specfreqs % vector of spectral frequencies 
%
%         clustinfo.erspdata    % (ncomps,ntimes,nfreqs) array of component ERSPs
%           clustinfo.ersptimes % vector of ERSP latencies
%           clustinfo.erspfreqs % vector of ERSP frequencies
%
%         clustinfo.itcdata     % (ncomps,ntimes,nfreqs) array of component ITCs
%           clustinfo.itctimes  % vector of ITC latencies
%           clustinfo.itc_freqs % vector of ITC frequencies
%
%  ********** SOON **********
%
%         clustinfo.topo        % (ncomps,65,65) array of component scalp map grids
%           clustinfo.topox     % abscissa values for columns of the scalp maps
%           clustinfo.topoy     % ordinate values for rows of the scalp maps
%
%         clustinfo.dipole      % array of component dipole information structs
%                               % with same format as EEG.dipfit.model
% Example:
%         % To obtain the ERPs for all Cluster-3 components from a STUDY
%         %
%         [STUDY clustinfo] = std_clustread(STUDY, ALLEEG, 3, 'erp');
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

function [STUDY, clustinfo, allinds] = std_readdata(STUDY, ALLEEG, varargin);

if nargin < 2
    help std_readdata;
    return;
end
if nargin < 3
    channel = {};
end

opt = finputcheck( varargin, { 'condition'  'cell'    []       {};
                             'subject'    'cell'    []       {};
                             'channels'   'cell'    []       {};
                             'clusters'   'integer' []       [];
                             'freqrange'  'real'    []       [];
                             'timerange'  'real'    []       [];
                             'group'      'cell'    []       {};
                             'statmode'   'string'  { 'subjects' 'individual' 'common' 'trials' }       'individual';
                             'subbaseline' 'string'  { 'on' 'off' }       'on';
                             'infotype'   'string'  { 'erp' 'spec' 'ersp' 'itc' } 'erp' }, 'std_readdata');
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

nc = length(STUDY.condition);
ng = length(STUDY.group);

% find channel indices
% --------------------
allinds = [];
if ~isempty(opt.channels)
    for c = 1:length(opt.channels)
        chanind = strmatch( lower(opt.channels{c}), lower({ STUDY.changrp.name }), 'exact');
        if isempty(chanind), error('Channel group not found'); end;
        allinds   = [ allinds chanind ];
    end;
else
    allinds = opt.clusters;
end;
    
% check if data is present
% ------------------------
for ind = 1:length(allinds)
    
    if ~isempty(opt.channels)
         tmpstruct = STUDY.changrp(allinds(ind));
         alldatasets = 1:length(STUDY.datasetinfo);
         allchanorcomp = -tmpstruct.chaninds;
    else tmpstruct = STUDY.cluster(allinds(ind));
         alldatasets = tmpstruct.sets; 
         allchanorcomp = repmat(tmpstruct.comps, [length(STUDY.condition) 1]);
    end;
    alldatasets   = alldatasets(:)';
    allchanorcomp = allchanorcomp(:)';
    
    dataread = 0;
    switch opt.infotype
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
                allerp  = cell( length(STUDY.condition), length(STUDY.group) );
                condgrp = cell( length(STUDY.condition), length(STUDY.group) );
                [ tmp alltimes ] = std_readerp( ALLEEG, 1, sign(allchanorcomp(1)), opt.timerange);
                for cond = 1:nc
                    for group = 1:ng
                        allerp{cond, group} = single(zeros(length(alltimes), length(STUDY.subject)));
                    end;
                end;

                % read the data and select channels
                % ---------------------------------
                fprintf('Reading ERP data:');
                for indtmp = 1:length(alldatasets)

                    index = alldatasets(indtmp);
                    condind = strmatch( STUDY.datasetinfo(index).condition, STUDY.condition ); if isempty(condind), condind = 1; end;
                    grpind  = strmatch( STUDY.datasetinfo(index).group    , STUDY.group     ); if isempty(grpind) , grpind  = 1; end;
                    subject = strmatch( STUDY.datasetinfo(index).subject  , STUDY.subject   ); 

                    % load data
                    % ---------
                    tmpind = allchanorcomp(indtmp); 
                    condgrp{condind, grpind}(subject) = tmpind;                    
                    if ~isnan(tmpind)
                        [ tmperp alltimes ] = std_readerp( ALLEEG, index, tmpind, opt.timerange);
                        allerp{condind, grpind}(:,subject) = tmperp(:);
                        fprintf('.');
                    else
                        allerp{condind, grpind}(:,subject) = NaN;
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
                allspec = cell( length(STUDY.condition), length(STUDY.group) );
                condgrp = cell( length(STUDY.condition), length(STUDY.group) );
                %[ tmpersp allfreqs alltimes tmpparams tmpspec] = std_readersp( ALLEEG, 1, -1, [], opt.freqrange);
                [ tmp allfreqs ] = std_readspec( ALLEEG, 1, sign(allchanorcomp(1)), opt.freqrange);
                for cond = 1:nc
                    for group = 1:ng
                        allspec{cond, group} = single(zeros(length(allfreqs), length(STUDY.subject)));
                    end;
                end;

                % read the data and select channels
                % ---------------------------------
                fprintf('Reading Spectrum data:');
                for indtmp = 1:length(alldatasets)

                    index = alldatasets(indtmp);
                    condind = strmatch( STUDY.datasetinfo(index).condition, STUDY.condition ); if isempty(condind), condind = 1; end;
                    grpind  = strmatch( STUDY.datasetinfo(index).group    , STUDY.group     ); if isempty(grpind) , grpind  = 1; end;
                    subject = strmatch( STUDY.datasetinfo(index).subject  , STUDY.subject   ); 
                    
                    % load data
                    % ---------
                    tmpind = allchanorcomp(indtmp);
                    condgrp{condind, grpind}(subject) = tmpind;                    
                    if ~isnan(tmpind)
                        %[ tmpersp allfreqs alltimes tmpparams tmpspec] = std_readersp( ALLEEG, index, tmpind, [], opt.freqrange);
                        %allspec{condind, grpind}(:,subject) = 10*log(tmpspec(:));
                        [ tmpspec allfreqs ] = std_readspec( ALLEEG, index, tmpind, opt.freqrange);
                        allspec{condind, grpind}(:,subject) = tmpspec(:);
                        fprintf('.');
                    else
                        allspec{condind, grpind}(:,subject) = NaN;
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
                ersp     = cell( length(STUDY.condition), length(STUDY.group) );
                erspbase = cell( length(STUDY.condition), length(STUDY.group) );
                erspinds = cell( length(STUDY.condition), length(STUDY.group) );
                condgrp = cell( length(STUDY.condition), length(STUDY.group) );
                [ tmp allfreqs alltimes ] = std_readersp( ALLEEG, 1, Inf*sign(allchanorcomp(1)), opt.timerange, opt.freqrange);
                for c = 1:nc
                    for g = 1:ng
                        ersp{c, g}     = single(zeros(length(alltimes), length(allfreqs), length(STUDY.subject)));
                        erspbase{c, g} = single(zeros(               1, length(allfreqs), length(STUDY.subject)));
                        if strcmpi(opt.statmode, 'trials')
                            ersp{c, g}  = single(zeros(length(alltimes), length(allfreqs), tottrials{c, g}));
                            count{c, g} = 1;
                        end;
                    end;
                end;
            
                % read the data and select channels
                % ---------------------------------
                fprintf('Reading all %s data:', upper(opt.infotype));
                for indtmp = 1:length(alldatasets)

                    index = alldatasets(indtmp);
                    condind = strmatch( STUDY.datasetinfo(index).condition, STUDY.condition ); if isempty(condind), condind = 1; end;
                    grpind  = strmatch( STUDY.datasetinfo(index).group    , STUDY.group     ); if isempty(grpind) , grpind  = 1; end;
                    subject = strmatch( STUDY.datasetinfo(index).subject  , STUDY.subject   ); 

                    % load data
                    % ---------
                    tmpind = allchanorcomp(indtmp);
                    condgrp{condind, grpind}(subject) = tmpind;                    
                    if strcmpi(opt.statmode, 'trials')
                        [ tmpersp allfreqs alltimes tmpparams] = std_readtimef( ALLEEG, index, tmpind, opt.timerange, opt.freqrange);
                        indices = [count{condind, grpind}:count{condind, grpind}+size(tmpersp,3)-1];
                        ersp {condind, grpind}(:,:,indices) = single(permute(tmpersp, [2 1 3]));
                        erspinds{condind, grpind}(1:2,subject) = [ count{condind, grpind} count{condind, grpind}+size(tmpersp,3)-1 ];
                        count{condind, grpind} = count{condind, grpind}+size(tmpersp,3);
                        if size(tmpersp,3) ~= ALLEEG(STUDY.datasetinfo(index).index).trials
                            error( sprintf('Wrong number of trials in datafile for dataset %d\n', STUDY.datasetinfo(index).index));
                        end;
                    elseif strcmpi(opt.infotype, 'itc')
                        [ tmpersp allfreqs alltimes tmpparams] = std_readitc( ALLEEG, index, tmpind, opt.timerange, opt.freqrange);
                        ersp{condind, grpind}(:,:,subject)     = single(abs(permute(tmpersp    , [2 1])));
                    else
                        [ tmpersp allfreqs alltimes tmpparams tmperspbase] = std_readersp( ALLEEG, index, tmpind, opt.timerange, opt.freqrange);
                        ersp{condind, grpind}(:,:,subject)     = single(permute(tmpersp    , [2 1]));
                        erspbase{condind, grpind}(:,:,subject) = single(10*log(permute(tmperspbase, [2 1])));
                    end;
                    %catch,
                    %    error(sprintf('Error loading time-freq file for dataset ''%s'': has it been computed?\n', ALLEEG.filename));
                    %end;
                    fprintf('.');

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
                    for g = 1:ng
                        for c = 1:nc
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
    end; % end switch
    if isempty(opt.channels) & exist('condgrp'), tmpstruct.compinds = condgrp; end;
    
    % copy results to structure
    % -------------------------
    fieldnames = { 'erpdata' 'erptimes' 'specdata' 'specfreqs' 'erspdata' 'erspbase' 'erspfreqs' 'ersptimes' ...
                   'itcfreqs' 'itctimes' 'itcdata' 'erspsubjinds' 'itcsubjinds' 'compinds' };
    for f = 1:length(fieldnames)
        if isfield(tmpstruct, fieldnames{f}), 
            tmpdata = getfield(tmpstruct, fieldnames{f});
            if ~isempty(opt.channels)
                 STUDY.changrp = setfield(STUDY.changrp, { allinds(ind) }, fieldnames{f}, tmpdata);
            else STUDY.cluster = setfield(STUDY.cluster, { allinds(ind) }, fieldnames{f}, tmpdata);
            end;
        end;
    end;
end;

% return structure
% ----------------
if ~isempty(opt.channels)
     clustinfo = STUDY.changrp(allinds);
else clustinfo = STUDY.cluster(allinds);
end;
