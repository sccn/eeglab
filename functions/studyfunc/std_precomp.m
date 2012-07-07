% std_precomp() - Precompute measures (ERP, spectrum, ERSP, ITC) for channels in a  
%                 study. If channels are interpolated before computing the measures,
%                 the updated EEG datasets are also saved to disk. Called by 
%                 pop_precomp(). Follow with pop_plotstudy(). See Example below.
% Usage:    
% >> [ALLEEG,STUDY] = std_precomp(STUDY, ALLEEG, chanorcomp, 'key', 'val', ...);
%
% Required inputs:
%   STUDY        - an EEGLAB STUDY set of loaded EEG structures
%   ALLEEG       - ALLEEG vector of one or more loaded EEG dataset structures
%   chanorcomp   - ['components'|'channels'| or channel cell array] The string 
%                  'components' forces the program to precompute all selected 
%                  measures for components. The string 'channels' forces the 
%                  program to compute all measures for all channels.
%                  A channel cell array containing channel labels will precompute
%                  the selected measures. Note that the name of the channel is
%                  not case-sensitive.
% Optional inputs:
%  'design'   - [integer] use specific study index design to compute measure.
%  'erp'      - ['on'|'off'] pre-compute ERPs for each dataset.
%  'spec'     - ['on'|'off'] pre-compute spectrum for each dataset.
%               Use 'specparams' to set spectrum parameters.
%  'ersp'     - ['on'|'off'] pre-compute ERSP for each dataset.
%               Use 'erspparams' to set time/frequency parameters.
%  'itc'      - ['on'|'off'] pre-compute ITC for each dataset.
%               Use 'erspparams' to set time/frequency parameters.
%  'scalp'    - ['on'|'off'] pre-compute scalp maps for components.
%  'allcomps' - ['on'|'off'] compute ERSP/ITC for all components ('off'
%               only use pre-selected components in the pop_study interface).
%  'erpparams'   - [cell array] Parameters for the std_erp function. See 
%                  std_erp for more information.
%  'specparams'  - [cell array] Parameters for the std_spec function. See 
%                  std_spec for more information.
%  'erspparams'  - [cell array] Optional arguments for the std_ersp function.
%  'erpimparams' - [cell array] Optional argument for std_erpimage. See
%                  std_erpimage for the list of arguments.
%  'recompute'   - ['on'|'off'] force recomputing ERP file even if it is 
%                  already on disk.
%  'rmicacomps'  - ['on'|'off'] remove ICA components pre-selected in each
%                  dataset (EEGLAB menu item, "Tools > Reject data using ICA 
%                  > Reject components by map). This option is ignored when 
%                  precomputing measures for ICA clusters. Default is 'off'.
%  'rmclust'     - [integer array] remove selected ICA component clusters.
%                  For example, ICA component clusters containing
%                  artifacts. This option is ignored when precomputing
%                  measures for ICA clusters.
%  'savetrials'  - ['on'|'off'] save single-trials ERSP. Requires a lot of disk
%                  space (dataset space on disk times 10) but allow for refined
%                  single-trial statistics.
% Outputs:
%   ALLEEG       - the input ALLEEG vector of EEG dataset structures, modified  
%                  by adding preprocessing data as pointers to Matlab files that 
%                  hold the pre-clustering component measures.
%   STUDY        - the input STUDY set with pre-clustering data added,
%                  for use by pop_clust()
%
% Example:
%   >> [STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, { 'cz' 'oz' }, 'interpo', ...
%               'on', 'erp', 'on', 'spec', 'on', 'ersp', 'on', 'erspparams', ...
%               { 'cycles' [ 3 0.5 ], 'alpha', 0.01, 'padratio' 1 });
%                          
%           % This prepares, channels 'cz' and 'oz' in the STUDY datasets.
%           % If a data channel is missing in one dataset, it will be
%           % interpolated (see eeg_interp()). The ERP, spectrum, ERSP, and 
%           % ITC for each dataset is then computed. 
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, 2006-

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2006, arno@sccn.ucsd.edu
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

function [ STUDY, ALLEEG ] = std_precomp(STUDY, ALLEEG, chanlist, varargin)
    
    if nargin < 2
        help std_precomp;
        return;
    end;
    
    if nargin == 2
        chanlist = 'channels'; % default to clustering the whole STUDY 
    end   
    Ncond = length(STUDY.condition);
    if Ncond == 0
        Ncond = 1;
    end

    g = finputcheck(varargin, { 'erp'         'string'  { 'on','off' }     'off';
                                'interp'      'string'  { 'on','off' }     'off';
                                'ersp'        'string'  { 'on','off' }     'off';
                                'recompute'   'string'  { 'on','off' }     'off';
                                'spec'        'string'  { 'on','off' }     'off';
                                'erpim'       'string'  { 'on','off' }     'off';
                                'scalp'       'string'  { 'on','off' }     'off';
                                'allcomps'    'string'  { 'on','off' }     'off';
                                'itc'         'string'  { 'on','off' }     'off';
                                'savetrials'  'string'  { 'on','off' }     'off';
                                'rmicacomps'  'string'  { 'on','off' }     'off';
                                'design'      'integer' []                 STUDY.currentdesign;
                                'rmclust'     'integer' []                 [];
                                'rmbase'      'integer' []                 []; % deprecated, for backward compatibility purposes, not documented
                                'specparams'        'cell'    {}                 {};
                                'erpparams'         'cell'    {}                 {};
                                'erpimparams'       'cell'    {}                 {};
                                'erspparams'        'cell'    {}                 {}}, 'std_precomp');
    if isstr(g), error(g); end;
    if ~isempty(g.rmbase), g.erpparams = { g.erpparams{:} 'rmbase' g.rmbase }; end;
    
    % union of all channel structures
    % -------------------------------
    computewhat = 'channels';
    if isstr(chanlist)
        if strcmpi(chanlist, 'channels')
            chanlist = [];
        else % components
            computewhat = 'components';
            if strcmpi(g.allcomps, 'on')
                chanlist = {};
                for index = 1:length(STUDY.datasetinfo)
                    chanlist = { chanlist{:} [1:size(ALLEEG(STUDY.datasetinfo(index).index).icaweights,1)] };
                end;
            else
                chanlist = { STUDY.datasetinfo.comps };
            end;
        end;
    end;
    if isempty(chanlist)
        alllocs = eeg_mergelocs(ALLEEG.chanlocs);
        chanlist = { alllocs.labels };
    elseif ~isnumeric(chanlist{1})
        alllocs = eeg_mergelocs(ALLEEG.chanlocs);
        [tmp c1 c2] = intersect( lower({ alllocs.labels }), lower(chanlist));
        [tmp c2] = sort(c2);
        alllocs = alllocs(c1(c2));
    end;
    
    % test if interp and reconstruct channel list
    % -------------------------------------------
    if strcmpi(computewhat, 'channels')
        if strcmpi(g.interp, 'on')
            STUDY.changrp = [];
            STUDY = std_changroup(STUDY, ALLEEG, chanlist, 'interp');
            g.interplocs = alllocs;
        else
            STUDY.changrp = [];
            STUDY = std_changroup(STUDY, ALLEEG, chanlist);
            g.interplocs = struct([]);
        end;
    end;
    
    % components or channels
    % ----------------------
    if strcmpi(computewhat, 'channels')
         curstruct = STUDY.changrp;
    else curstruct = STUDY.cluster;
    end;
    
    % compute ERPs
    % ------------
    if strcmpi(g.erp, 'on')
        % check dataset consistency
        % -------------------------
        allPnts = [ALLEEG([STUDY.design(g.design).cell.dataset]).pnts];
        if iscell(allPnts), allPnts = [ allPnts{:} ]; end;
        if length(unique(allPnts)) > 1
            error([ 'Cannot compute ERPs because datasets' 10 'do not have the same number of data points' ])
        end;
        
        for index = 1:length(STUDY.design(g.design).cell)
            desset = STUDY.design(g.design).cell(index);
            addopts = { 'savetrials', g.savetrials, 'recompute', g.recompute, 'fileout', desset.filebase, 'trialindices', desset.trials };
            if strcmpi(computewhat, 'channels')
                [tmpchanlist opts] = getchansandopts(STUDY, ALLEEG, chanlist, desset.dataset, g);
                std_erp(ALLEEG(desset.dataset), 'channels', tmpchanlist, opts{:}, addopts{:}, g.erpparams{:});
            else
                if length(desset.dataset)>1 && ~isequal(chanlist{desset.dataset})
                    error(['ICA decompositions must be identical if' 10 'several datasets are concatenated to build' 10 'the design, abording' ]);
                end;
                std_erp(ALLEEG(desset.dataset), 'components', chanlist{desset.dataset(1)}, addopts{:}, g.erpparams{:});
            end;
        end;
        if isfield(curstruct, 'erpdata')
            curstruct = rmfield(curstruct, 'erpdata');
            curstruct = rmfield(curstruct, 'erptimes');
        end;
    end;
    
    % compute spectrum
    % ----------------
    if strcmpi(g.spec, 'on')
        % check dataset consistency
        % -------------------------
        for index = 1:length(STUDY.design(g.design).cell)
            desset = STUDY.design(g.design).cell(index);
            addopts = { 'savetrials', g.savetrials, 'recompute', g.recompute, 'fileout', desset.filebase, 'trialindices', desset.trials };
            if strcmpi(computewhat, 'channels')
                [tmpchanlist opts] = getchansandopts(STUDY, ALLEEG, chanlist, desset.dataset, g);
                std_spec(ALLEEG(desset.dataset), 'channels', tmpchanlist, opts{:}, addopts{:}, g.specparams{:});
            else
                if length(desset.dataset)>1 && ~isequal(chanlist{desset.dataset})
                    error(['ICA decompositions must be identical if' 10 'several datasets are concatenated to build' 10 'the design, abording' ]);
                end;
                std_spec(ALLEEG(desset.dataset), 'components', chanlist{desset.dataset(1)}, addopts{:}, g.specparams{:});
            end;
        end;
        if isfield(curstruct, 'specdata')
            curstruct = rmfield(curstruct, 'specdata');
            curstruct = rmfield(curstruct, 'specfreqs');
        end;
    end;

    % compute spectrum
    % ----------------
    if strcmpi(g.erpim, 'on')
        % check dataset consistency
        % -------------------------
        allPnts = [ALLEEG([STUDY.design(g.design).cell.dataset]).pnts];
        if iscell(allPnts), allPnts = [ allPnts{:} ]; end;
        if length(unique(allPnts)) > 1
            error([ 'Cannot compute ERSPs/ITCs because datasets' 10 'do not have the same number of data points' ])
        end;
        
        % check consistency with parameters on disk
        % -----------------------------------------
        guimode = 'guion';
        tmpparams = {};
        tmpparams = g.erpimparams;
        if strcmpi(g.recompute, 'off')
            for index = 1:length(STUDY.design(g.design).cell)
                desset = STUDY.design(g.design).cell(index);
                if strcmpi(computewhat, 'channels')
                    filename = [ desset.filebase '.daterpim'];
                else filename = [ desset.filebase '.icaerpim'];
                end;
                [guimode, g.erpimparams] = std_filecheck(filename, g.erpimparams, guimode, { 'fileout' 'recompute', 'channels', 'components', 'trialindices'});
                if strcmpi(guimode, 'cancel'), return; end;
            end;
            if strcmpi(guimode, 'usedisk') || strcmpi(guimode, 'same'), g.recompute = 'off';
            else                                                        g.recompute = 'on';
            end;
            if ~isempty(g.erpimparams) && isstruct(g.erpimparams)
                tmpparams      = fieldnames(g.erpimparams); tmpparams = tmpparams';
                tmpparams(2,:) = struct2cell(g.erpimparams);
            end;
        end;
        
        % set parameters in ERPimage parameters
        % -------------------------------------
        STUDY = pop_erpimparams(STUDY, tmpparams{:}); % a little trashy as the function pop_erpimparams does not check the fields
        
        % compute ERPimages
        % -----------------
        for index = 1:length(STUDY.design(g.design).cell)
            desset = STUDY.design(g.design).cell(index);
            addopts = { 'recompute', g.recompute, 'fileout', desset.filebase, 'trialindices', desset.trials };
            if strcmpi(computewhat, 'channels')
                [tmpchanlist opts] = getchansandopts(STUDY, ALLEEG, chanlist, desset.dataset, g);
                std_erpimage(ALLEEG(desset.dataset), 'channels', tmpchanlist, opts{:}, addopts{:}, tmpparams{:});
            else
                if length(desset.dataset)>1 && ~isequal(chanlist{desset.dataset})
                    error(['ICA decompositions must be identical if' 10 'several datasets are concatenated to build' 10 'the design, abording' ]);
                end;
                std_erpimage(ALLEEG(desset.dataset), 'components', chanlist{desset.dataset(1)}, addopts{:}, tmpparams{:});
            end;
        end;
        if isfield(curstruct, 'erpimdata')
            curstruct = rmfield(curstruct, 'erpimdata');
            curstruct = rmfield(curstruct, 'erpimtimes');
            curstruct = rmfield(curstruct, 'erpimtrials');
            curstruct = rmfield(curstruct, 'erpimevents');
        end;
    end;
    
    % compute component scalp maps
    % ----------------------------
    if strcmpi(g.scalp, 'on') && ~strcmpi(computewhat, 'channels')
        for index = 1:length(STUDY.datasetinfo)
            
            % find duplicate
            % --------------
            found = [];
            ind1 = STUDY.datasetinfo(index).index;
            inds = strmatch(STUDY.datasetinfo(index).subject, { STUDY.datasetinfo(1:index-1).subject });
            for index2 = 1:length(inds)
                ind2 = STUDY.datasetinfo(inds(index2)).index;
                if isequal(ALLEEG(ind1).icawinv, ALLEEG(ind2).icawinv)
                    found = ind2;
                end;
            end;
            
            % make link if duplicate
            % ----------------------
            fprintf('Computing/checking topo file for dataset %d\n', ind1);
            if ~isempty(found)
                tmpfile1 = fullfile( ALLEEG(index).filepath, [ ALLEEG(index).filename(1:end-3) 'icatopo' ]); 
                tmp.file = fullfile( ALLEEG(found).filepath, [ ALLEEG(found).filename(1:end-3) 'icatopo' ]); 
                std_savedat(tmpfile1, tmp);
            else
                std_topo(ALLEEG(index), chanlist{index}, 'none', 'recompute', g.recompute);
            end;
        end;
        if isfield(curstruct, 'topo')
            curstruct = rmfield(curstruct, 'topo');
            curstruct = rmfield(curstruct, 'topox');
            curstruct = rmfield(curstruct, 'topoy');
            curstruct = rmfield(curstruct, 'topoall');
            curstruct = rmfield(curstruct, 'topopol');
        end;
    end;
    
    % compute ERSP and ITC
    % --------------------
    if strcmpi(g.ersp, 'on') || strcmpi(g.itc, 'on')
        % check dataset consistency
        % -------------------------
        allPnts = [ALLEEG([STUDY.design(g.design).cell.dataset]).pnts];
        if iscell(allPnts), allPnts = [ allPnts{:} ]; end;
        if length(unique(allPnts)) > 1
            error([ 'Cannot compute ERSPs/ITCs because datasets' 10 'do not have the same number of data points' ])
        end;
        
        if strcmpi(g.ersp, 'on') & strcmpi(g.itc, 'on'), type = 'both';
        elseif strcmpi(g.ersp, 'on')                   , type = 'ersp';
        else                                             type = 'itc';
        end;
        
        % check for existing files
        % ------------------------
        guimode = 'guion';
        [ tmpX tmpt tmpf g.erspparams ] = std_ersp(ALLEEG(1), 'channels', 1, 'type', type, 'recompute', 'on', 'getparams', 'on', 'savetrials', g.savetrials, g.erspparams{:});
        if strcmpi(g.recompute, 'off')
            for index = 1:length(STUDY.design(g.design).cell)
                desset = STUDY.design(g.design).cell(index);
                if strcmpi(computewhat, 'channels')
                     filename = [ desset.filebase '.datersp'];
                else filename = [ desset.filebase '.icaersp'];
                end;
                [guimode, g.erspparams] = std_filecheck(filename, g.erspparams, guimode, { 'plotitc' 'plotersp' 'plotphase' });
                if strcmpi(guimode, 'cancel'), return; end;
            end;
            if strcmpi(guimode, 'usedisk') || strcmpi(guimode, 'same'), g.recompute = 'off'; 
            else                                                        g.recompute = 'on'; 
            end;
        end;
        
        % check for existing files
        % ------------------------
        if isempty(g.erspparams), 
            tmpparams = {};
        elseif iscell(g.erspparams), 
            tmpparams = g.erspparams; 
        else
            tmpparams      = fieldnames(g.erspparams); tmpparams = tmpparams';
            tmpparams(2,:) = struct2cell(g.erspparams);
        end;
        tmpparams = { tmpparams{:} 'recompute' g.recompute };
        for index = 1:length(STUDY.design(g.design).cell)
            desset = STUDY.design(g.design).cell(index);
            if strcmpi(computewhat, 'channels')
                [tmpchanlist opts] = getchansandopts(STUDY, ALLEEG, chanlist, desset.dataset, g);
                std_ersp(ALLEEG(desset.dataset), 'channels', tmpchanlist, 'type', type, 'fileout', desset.filebase, 'trialindices', desset.trials, opts{:}, tmpparams{:});
            else
                if length(desset.dataset)>1 && ~isequal(chanlist{desset.dataset})
                    error(['ICA decompositions must be identical if' 10 'several datasets are concatenated to build' 10 'the design, abording' ]);
                end;
                std_ersp(ALLEEG(desset.dataset), 'components', chanlist{desset.dataset(1)}, 'type', type, 'fileout', desset.filebase, 'trialindices', desset.trials, tmpparams{:});
            end;
        end;
        if isfield(curstruct, 'erspdata')
            curstruct = rmfield(curstruct, 'erspdata');
            curstruct = rmfield(curstruct, 'ersptimes');
            curstruct = rmfield(curstruct, 'erspfreqs');
        end;
        if isfield(curstruct, 'itcdata')
            curstruct = rmfield(curstruct, 'itcdata');
            curstruct = rmfield(curstruct, 'itctimes');
            curstruct = rmfield(curstruct, 'itcfreqs');
        end;
    end;

    % components or channels
    % ----------------------
    if strcmpi(computewhat, 'channels')
         STUDY.changrp = curstruct;
    else STUDY.cluster = curstruct;
    end;
    
    return;
        
    % find components in cluster for specific dataset
    % -----------------------------------------------
    function rmcomps = getclustcomps(STUDY, rmclust, settmpind);    
    
        rmcomps   = cell(1,length(settmpind));
        for idat = 1:length(settmpind) % scan dataset for which to find component clusters
            currentdat = settmpind(idat);
            for rmi = 1:length(rmclust) % scan clusters
                compinds = STUDY.cluster(rmclust(rmi)).allinds;
                
                for ind = 1:length(STUDY.cluster(rmclust(rmi)).setinds(:)) % scan datasets of cluster
                    setinds  = [STUDY.design(STUDY.currentdesign).cell(STUDY.cluster(rmclust(rmi)).setinds{ind}).dataset ];
                    indmatch = find(setinds == currentdat);
                    if ~isempty(indmatch)
                        rmcomps{idat} = union( rmcomps{idat}, compinds{ind}(indmatch));
                    end;
                end;
            end;
        end;
        
    % make option array and channel list (which depend on interp) for any type of measure
    % ----------------------------------------------------------------------
    function [tmpchanlist, opts] = getchansandopts(STUDY, ALLEEG, chanlist, idat, g);
        
        opts = { };
        if ~isempty(g.rmclust)
            opts = { opts{:} 'rmcomps' getclustcomps(STUDY, g.rmclust, idat) };                
        elseif strcmpi(g.rmicacomps, 'on')
            for ind = 1:length(idat)
                rmcomps{ind} = find(ALLEEG(idat(1)).reject.gcompreject);
            end;
            opts = { opts{:} 'rmcomps' rmcomps };
        end;
        if strcmpi(g.interp, 'on')
            tmpchanlist = chanlist;
            allocs = eeg_mergelocs(ALLEEG.chanlocs);
            [tmp1 tmp2 neworder] = intersect( {allocs.labels}, chanlist);
            [tmp1 ordertmp2] = sort(tmp2);
            neworder = neworder(ordertmp2);
            opts = { opts{:} 'interp' allocs(neworder) };
        else
            newchanlist = [];
            tmpchanlocs = ALLEEG(idat(1)).chanlocs;
            chanlocs = { tmpchanlocs.labels };
            for i=1:length(chanlist)
                newchanlist = [ newchanlist strmatch(chanlist{i}, chanlocs, 'exact') ];
            end;
            tmpchanlocs =  ALLEEG(idat(1)).chanlocs;
            tmpchanlist = { tmpchanlocs(newchanlist).labels };
        end;
        
