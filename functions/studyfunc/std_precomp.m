% std_precomp() - Precompute measures (ERP, spectrum, ERSP, ITC) for channels in a study. 
%                 If channels are interpolated before computing the measures, the updated 
%                 EEG datasets are also saved to disk. Called by pop_precomp(). Follow with 
%                 pop_plotstudy(). See Example below.
% Usage:    
% >> [ALLEEG,STUDY] = std_precomp(ALLEEG, STUDY, chanlist, 'key', 'val', ...);
%
% Required inputs:
%   ALLEEG       - ALLEEG vector of one or more loaded EEG dataset structures
%   STUDY        - an EEGLAB STUDY set of loaded EEG structures
%   chanlist     - [cell array] Channel names for which to precompute the
%                  selected measures. Note that the name of the channel is
%                  not case-sensitive.
% Optional inputs:
%  'erp'     - ['on'|'off'] pre-compute ERPs for each dataset.
%  'spec'    - ['on'|'off'] pre-compute spectrum for each dataset.
%              Use 'specparams' to set spectrum parameters.
%  'ersp'    - ['on'|'off'] pre-compute ERSP for each dataset.
%              Use 'erspparams' to set time/frequency parameters.
%  'itc'     - ['on'|'off'] pre-compute ITC for each dataset.
%              Use 'erspparams' to set time/frequency parameters.
%  'specparams' - [cell array] Parameters for the spectopo function are given as 
%              optional arguments. Note that it is advised to compute spectrum 
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
        chanlist = []; % default to clustering the whole STUDY 
    end    
    Ncond = length(STUDY.condition);
    if Ncond == 0
        Ncond = 1;
    end

    g = finputcheck(varargin, { 'erp'         'string'  { 'on' 'off' }     'off';
                                'interpolate' 'string'  { 'on' 'off' }     'off';
                                'ersp'        'string'  { 'on' 'off' }     'off';
                                'spec'        'string'  { 'on' 'off' }     'off';
                                'itc'         'string'  { 'on' 'off' }     'off';
                                'specparams'        'cell'    {}                 {};
                                'erspparams'        'cell'    {}                 {}}, 'std_precomp');
    if isstr(g), error(g); end;
    
    % union of all channel structures
    % -------------------------------
    if isempty(chanlist)
        alllocs = ALLEEG(STUDY.datasetinfo(1).index).chanlocs;
        alllabs = { alllocs.labels };
        for index = 2:length(STUDY.datasetinfo)
           tmplocs = ALLEEG(STUDY.datasetinfo(index).index).chanlocs;
           alllocs = eeg_mergechan(alllocs, tmplocs);
        end;
        chanlist = { alllocs.labels };
    end;
    
    % Interpolate all datasets first
    % ------------------------------
    if strcmpi(g.interpolate, 'on')
        fprintf('Interpolation of data channels\n');
        fprintf('------------------------------\n');
        [ STUDY, ALLEEG ] = std_interp(STUDY, ALLEEG, chanlist);
    end;
    
    % compute ERPs
    % ------------
    if strcmpi(g.erp, 'on')
        for index = 1:length(STUDY.datasetinfo)
            std_erp(ALLEEG(STUDY.datasetinfo(index).index), 'channels', chanlist);
        end;
    end;

    % compute spectrum
    % ----------------
    if strcmpi(g.spec, 'on')
        for index = 1:length(STUDY.datasetinfo)
            std_spec(ALLEEG(STUDY.datasetinfo(index).index), 'channels', chanlist, g.specparams{:});
        end;
    end;

    % compute ERSP
    % ------------
    if strcmpi(g.ersp, 'on') | strcmpi(g.itc, 'on')
        if strcmpi(g.ersp, 'on') & strcmpi(g.itc, 'on'), type = 'both';
        elseif strcmpi(g.ersp, 'on')                   , type = 'ersp';
        else                                             type = 'itc';
        end;
        
        % check for existing files
        % ------------------------
        guimode = 'guion';
        for index = 1:length(STUDY.datasetinfo)
            
            filename = fullfile( ALLEEG(index).filepath,[ ALLEEG(index).filename(1:end-3) 'datersp']);
            [guimode, g.erspparams] = std_filecheck(filename, g.erspparams, guimode, { 'plotitc' 'plotersp' 'plotphase' });
            if strcmpi(guimode, 'cancel'), return; end;
            
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
            std_ersp(ALLEEG(STUDY.datasetinfo(index).index), 'channels', chanlist, 'type', type, tmpparams{:});
        end;
    end;
    
    STUDY = std_changroup(STUDY, ALLEEG);
    return;
