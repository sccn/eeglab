% pop_limo() - prepare and convert EEGLAB data and structure to be
%              processed by LIMO.
%
% Usage:
%   >> STUDY = pop_limo( STUDY, ALLEEG, 'key', val)
%
% Inputs:
%   STUDY        - an EEGLAB STUDY set of loaded EEG structures
%   ALLEEG       - ALLEEG vector of one or more loaded EEG dataset structures
%
% Optional inputs: same as std_limo()
%
% Outputs:
%   STUDY       - an EEGLAB STUDY set of loaded EEG structures
%
% Author: Arnaud Delorme, SCCN, UCSD, 2015-
%
% See also: std_limo()

% Copyright (C) Arnaud Delorme
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

function [STUDY,com] = pop_limo(STUDY, ALLEEG, varargin)

if nargin < 2
    help pop_limo;
    return;
end;
com = '';
if nargin < 3
    dataMeasures = { 'ERP' 'Spectrum' };
    fileMeasures = { 'daterp' 'datspec' };
    methods      = { 'OLS' 'WLS' };

    uilist = { ...
        {'style' 'text'       'string' [ 'LInear MOdeling of EEG data of design' int2str(STUDY.currentdesign) ] 'fontweight' 'bold' 'fontsize', 12} ...
        {'style' 'text'       'string' [ '(Use STUDY design interface to switch to a different design)' ] } {} ...
        {'style' 'text'       'string' 'Input data to use for the GLM' } ...
        {'style' 'popupmenu'  'string' dataMeasures 'tag' 'measure' } ...
        {'style' 'text'       'string' 'Optimization method' } ...
        {'style' 'popupmenu'  'string' methods 'tag' 'method' } ...
        {'style' 'checkbox'   'string' 'Erase previous model for this design and measure' 'tag' 'erase' 'value' true } ...
        };
    
    cline = [1.1 1.1];
    geometry = { 1 1 1 cline cline 1 };
    geomvert = [ 1 1 1 1     1     1 ];
        
    [out_param userdat tmp res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'geomvert', geomvert, ...
                                            'title', 'LInear MOdeling of EEG data -- pop_limo()');
    if isempty(res), return; end;
    
    options = { 'method' methods{res.method} 'measure' fileMeasures{res.measure} 'erase' fastif(res.erase, 'on', 'off') };
else
    options = varargin;
end;
if strcmp(fastif(res.measure, 'on', 'off'), 'on')
    if exist([STUDY.filepath filesep 'limo_batch_report'],'dir') == 7
        rmdir([STUDY.filepath filesep 'limo_batch_report'],'s');
    end
end
[STUDY tmp] = std_limo(STUDY, ALLEEG, options{:});
com = sprintf('pop_limo(STUDY, ALLEEG, %s);', vararg2str(options));
