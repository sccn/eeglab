% std_precomp() - Precompute measures (ERP, spectrum, ERSP, ITC) for channels in a study. 
%                 If channels are interpolated before computing the measures, the updated 
%                 EEG datasets are also saved to disk. Called by pop_precomp(). Follow with 
%                 pop_plotstudy(). See Example below.
% Usage:    
% >> [ALLEEG,STUDY] = std_precomp(ALLEEG, STUDY, chanorcomp, 'key', 'val', ...);
%
% Required inputs:
%   ALLEEG       - ALLEEG vector of one or more loaded EEG dataset structures
%   STUDY        - an EEGLAB STUDY set of loaded EEG structures
%   chanorcomp   - ['components'|'channels'| or channel cell array] The string 
%                  'components' forces the program to precompute all selected 
%                  measures for components. The string 'channels' forces the 
%                  program to compute all measures for all channels.
%                  A channel cell array containing channel labels will precompute
%                  the selected measures. Note that the name of the channel is
%                  not case-sensitive.
% Optional inputs:
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
%  'specparams' - [cell array] Parameters for the spectopo function are given as 
%              optional arguments:
%                   'freqrange' = [min max] frequency range to calculate. Changes x-axis limits {default: 
%                                1 Hz for the min and Nyquist (srate/2) for the max. If specified 
%                                power distribution maps are plotted, the highest mapped frequency 
%                                determines the max freq}. Note that it is better here to compute
%                                spectrum over a wide range of frequencies (it will then
%                                be possible to select another subrange for plotting).
%                   'freqfac'  = [integer] ntimes to oversample -> frequency resolution {default: 2}
%                   'nfft'     = [integer] length to zero-pad data to. Overwrites 'freqfac' above.
%                   'winsize'  = [integer] window size in data points {default: from data}
%                   'overlap'  = [integer] window overlap in data points {default: 0}
%                   'percent'  = [float 0 to 100] percent of the data to sample for computing the 
%                                spectra. Values < 100 speed up the computation. {default: 100}.
%                   'mapnorm'  = [float vector] If 'data' contain the activity of an independant 
%                                component, this parameter should contain its scalp map. In this case
%                                the spectrum amplitude will be scaled to component RMS scalp power.
%                                Useful for comparing component strengths {default: none}
%                   'rmdc'     = ['on'|'off'] 'on' -> remove DC {default: 'off'}  
%
%              Note that it is advised to compute spectrum 
%              over all frequencies since plotting function can always reduce
%              the range of plotted frequencies.
%  'erspparams' - [cell array] Optional arguments are 'cycles', 'freqrange',
%              'padratio', 'winsize', 'alpha' (see newtimef()). Note that it 
%              is adivised to select the largest frequency range and time window
%              as plotting function are capable of plotting subranges of
%              these. An important optional parameter that is
%                    'savetrials' = ['on'|'off'] save single-trials ERSP.
%                                   Requires a lot of disk space (dataset
%                                   space on disk times 10) but allow for
%                                   refined single-trial statistics.
%  'recompute' - ['on'|'off'] force recomputing ERP file even if it is 
%                already on disk.
%
% Outputs:
%   ALLEEG       - the input ALLEEG vector of EEG dataset structures, modified by adding preprocessing 
%                  data as pointers to Matlab files that hold the pre-clustering component measures.
%   STUDY        - the input STUDY set with pre-clustering data added, for use by pop_clust() 
%
% Example:
%   >> [ALLEEG STUDY] = std_precomp(ALLEEG, STUDY, { 'cz' 'oz' }, 'interpolate', 'on', 'erp', 'on', ...
%          'spec', 'on', 'ersp', 'on', 'erspparams', { 'cycles' [ 3 0.5 ], 'alpha', 0.01, 'padratio' 1 });
%                          
%           % This prepares, channels 'cz' and 'oz' in the STUDY datasets.
%           % If a data channel is missing in one dataset, it will be
%           % interpolated (see eeg_interp()). The ERP, spectrum, ERSP, and 
%           % ITC for each dataset is then computed. 
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, 2006-

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.21  2008/02/15 16:51:36  arno
% simplify code for merging channel location files
%
% Revision 1.20  2007/12/09 00:40:15  arno
% recompute for topo
%
% Revision 1.19  2007/11/22 23:34:55  arno
% header
%
% Revision 1.18  2007/11/22 23:34:03  arno
% header
%
% Revision 1.17  2007/11/21 16:40:53  arno
% help msg
%
% Revision 1.16  2007/10/25 00:59:39  nima
% spectopo parameters described in help message.
%
% Revision 1.15  2007/09/11 10:51:16  arno
% precompute measures for components
%
% Revision 1.14  2007/04/06 22:09:44  arno
% recompute tag
%
% Revision 1.13  2007/04/05 23:17:55  arno
% guimode
%
% Revision 1.12  2007/04/05 23:13:39  arno
% *** empty log message ***
%
% Revision 1.11  2007/02/28 12:05:14  arno
% option to force recomputation
%
% Revision 1.6  2007/01/29 10:50:27  arno
% fix ERSP options
%
% Revision 1.4  2006/11/14 04:12:53  arno
% [Asame
%
% Revision 1.3  2006/11/14 03:59:25  arno
% debug ERSP check
%
% Revision 1.2  2006/11/14 03:53:18  arno
% Now checking file on disk
%
% Revision 1.1  2006/09/12 18:43:54  arno
% Initial revision
%

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

    g = finputcheck(varargin, { 'erp'         'string'  { 'on' 'off' }     'off';
                                'interp'      'string'  { 'on' 'off' }     'off';
                                'ersp'        'string'  { 'on' 'off' }     'off';
                                'recompute'   'string'  { 'on' 'off' }     'off';
                                'spec'        'string'  { 'on' 'off' }     'off';
                                'scalp'       'string'  { 'on' 'off' }     'off';
                                'allcomps'    'string'  { 'on' 'off' }     'off';
                                'itc'         'string'  { 'on' 'off' }     'off';
                                'rmicacomps'  'string'  { 'on' 'off' }     'off';
                                'rmclust'     'integer' []                 [];
                                'rmbase'      'integer' []                 [];
                                'specparams'        'cell'    {}                 {};
                                'erspparams'        'cell'    {}                 {}}, 'std_precomp');
    if isstr(g), error(g); end;
    
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
        alllocs = eeg_mergelocs(ALLEEG(:).chanlocs);
        chanlist = { alllocs.labels };
    elseif ~isnumeric(chanlist{1})
        alllocs = eeg_mergelocs(ALLEEG(:).chanlocs);
        [tmp c1 c2] = intersect( lower({ alllocs.labels }), lower(chanlist));
        [tmp c2] = sort(c2);
        alllocs = alllocs(c1(c2));
    end;
    
    % test if interp and reconstruct channel list
    % -------------------------------------------
    if strcmpi(g.interp, 'on')
        STUDY = std_changroup(STUDY, ALLEEG, chanlist, 'interp');
        g.interplocs = alllocs;
    else
        STUDY = std_changroup(STUDY, ALLEEG, chanlist);
        g.interplocs = [];
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
        for index = 1:length(STUDY.datasetinfo)
            if strcmpi(computewhat, 'channels')
                [tmpchanlist opts] = getchansandopts(STUDY, ALLEEG, chanlist, index, g);
                std_erp(ALLEEG(STUDY.datasetinfo(index).index), 'channels', tmpchanlist, 'rmbase', g.rmbase, opts{:});
            else
                std_erp(ALLEEG(STUDY.datasetinfo(index).index), 'components', chanlist{index}, 'rmbase', g.rmbase, 'recompute', g.recompute);
            end;
        end;
        if isfield(curstruct, 'erpdata')
            curstruct = rmfield(curstruct, 'erpdata');
            curstruct = rmfield(curstruct, 'erptimes');
        end;
    end;

    % compute component scalp maps
    % ----------------------------
    if strcmpi(g.scalp, 'on')
        for index = 1:length(STUDY.datasetinfo)
            
            % find duplicate
            % --------------
            found = [];
            ind1 = STUDY.datasetinfo(index).index;
            inds = strmatch(STUDY.datasetinfo(index).subject, { STUDY.datasetinfo(1:index-1).subject });
            for index2 = inds'
                ind2 = STUDY.datasetinfo(index2).index;
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
    
    % compute spectrum
    % ----------------
    if strcmpi(g.spec, 'on')
        for index = 1:length(STUDY.datasetinfo)
            if strcmpi(computewhat, 'channels')
                [tmpchanlist opts] = getchansandopts(STUDY, ALLEEG, chanlist, index, g);
                std_spec(ALLEEG(STUDY.datasetinfo(index).index), 'channels', tmpchanlist, opts{:}, g.specparams{:});
            else
                std_spec(ALLEEG(STUDY.datasetinfo(index).index), 'components', chanlist{index}, 'recompute', g.recompute);
            end;
        end;
        if isfield(curstruct, 'specdata')
            curstruct = rmfield(curstruct, 'specdata');
            curstruct = rmfield(curstruct, 'specfreqs');
        end;
    end;

    % compute ERSP and ITC
    % --------------------
    if strcmpi(g.ersp, 'on') | strcmpi(g.itc, 'on')
        if strcmpi(g.ersp, 'on') & strcmpi(g.itc, 'on'), type = 'both';
        elseif strcmpi(g.ersp, 'on')                   , type = 'ersp';
        else                                             type = 'itc';
        end;
        
        % check for existing files
        % ------------------------
        guimode = 'guion';
        [ tmpX tmpt tmpf g.erspparams ] = std_ersp(ALLEEG(1), 'channels', 1, 'type', type, 'recompute', 'on', 'getparams', 'on', g.erspparams{:});
        if strcmpi(g.recompute, 'off')
            for index = 1:length(STUDY.datasetinfo)
            
                if strcmpi(computewhat, 'channels')
                    filename = fullfile( ALLEEG(index).filepath,[ ALLEEG(index).filename(1:end-3) 'datersp']);
                else
                    filename = fullfile( ALLEEG(index).filepath,[ ALLEEG(index).filename(1:end-3) 'icaersp']);
                end;
                [guimode, g.erspparams] = std_filecheck(filename, g.erspparams, guimode, { 'plotitc' 'plotersp' 'plotphase' });
                if strcmpi(guimode, 'cancel'), return; end;

            end;
            if strcmpi(guimode, 'usedisk') | strcmpi(guimode, 'same'), g.recompute = 'off'; 
            else                                                       g.recompute = 'on'; 
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
        for index = 1:length(STUDY.datasetinfo)
            if strcmpi(computewhat, 'channels')
                [tmpchanlist opts] = getchansandopts(STUDY, ALLEEG, chanlist, index, g);
                std_ersp(ALLEEG(STUDY.datasetinfo(index).index), 'channels', tmpchanlist, 'type', type, opts{:}, tmpparams{:});
            else
                std_ersp(ALLEEG(STUDY.datasetinfo(index).index), 'components', chanlist{index}, 'type', type, 'recompute', g.recompute, tmpparams{:});
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

    % get channel indices from changrp structure
    % ------------------------------------------    
    function chaninds = getchannelindices(changrp, datasetind);    
    
    chaninds = [];
    for index = 1:length(changrp)    
        tmpind = find([ changrp(index).setinds{:} ] == datasetind);
        if ~isempty(tmpind)
            chaninds = [ chaninds datasetind ];
        end;
    end;
        
    % find components in cluster for specific dataset
    % -----------------------------------------------
    function rmcomps = getclustcomps(STUDY, rmclust, settmpind);    
    
        rmcomps   = [ ];
        for rmi = 1:length(rmclust)
            findind   = find(settmpind == STUDY.cluster(rmclust(rmi)).setinds{c,g});
            rmcomps   = [ rmcomps STUDY.cluster(rmclust(rmi)).allinds{c,g}(findind) ];
        end;

    % make option array and channel list (which depend on interp) for any type of measure
    % ----------------------------------------------------------------------
    function [tmpchanlist opts] = getchansandopts(STUDY, ALLEEG, chanlist, index, g);
        
        idat = STUDY.datasetinfo(index).index;
        opts = { 'recompute' g.recompute };
        if ~isempty(g.rmclust)
            opts = { opts{:} 'rmcomps' getclustcomps(STUDY, g.rmclust, idat) };                
        elseif strcmpi(g.rmicacomps, 'on')
            opts = { opts{:} 'rmcomps' find(ALLEEG(idat).reject.gcompreject) };
        end;
        if ~isempty(g.interplocs)
            alllocs = eeg_mergelocs(ALLEEG(:).chanlocs);
            tmpchanlist = chanlist;
            opts = { opts{:} 'interp' g.interplocs };
        else
            tmpchanlist = getchannelindices(STUDY.changrp, STUDY.datasetinfo(index).index);
            tmpchanlist = { ALLEEG(STUDY.datasetinfo(index).index).chanlocs(tmpchanlist).labels };
        end;
        