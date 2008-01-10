% pop_erpimparams() - Set plotting and statistics parameters for cluster ERP 
%                   plotting
% Usage:    
%   >> STUDY = pop_erpimparams(STUDY, 'key', 'val');   
%
% Inputs:
%   STUDY        - EEGLAB STUDY set
%
% Statistics options:
%   'groupstats  - ['on'|'off'] Compute (or not) statistics across groups.
%                  {default: 'off'}
%   'condstats'  - ['on'|'off'] Compute (or not) statistics across groups.
%                  {default: 'off'}
%   'topotime'   - [real] Plot ERP scalp maps at one specific latency (ms).
%                   A latency range [min max] may also be defined (the 
%                   ERP is then averaged over the interval) {default: []}
%   'statistics' - ['param'|'perm'] Type of statistics to use: 'param' for
%                  parametric and 'perm' for permutation-based statistics. 
%                  {default: 'param'}
%   'naccu'      - [integer] Number of surrogate averages to accumulate for
%                  permutation statistics. For example, to test whether 
%                  p<0.01, use >=200. For p<0.001, use 'naccu' >=2000. 
%                  If a threshold (not NaN) is set below, and 'naccu' is 
%                  too low, it will be automatically reset. (This option 
%                  is now available only from the command line).
%   'threshold'  - [NaN|float<<1] Significance probability threshold. 
%                  NaN -> plot the p-values themselves on a different axis. 
%                  When possible, the significant time regions are indicated 
%                  below the data.
% See also: std_erpimage()
%
% Authors: Arnaud Delorme, CERCO, CNRS, 2008-

% Copyright (C) Arnaud Delorme, 2008
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
% Revision 1.2  2008/01/10 20:17:23  arno
% fixing sortvar -> sorttype
%
% Revision 1.1  2008/01/10 20:13:36  arno
% Initial revision
%

function [ STUDY, com ] = pop_erpimparams(STUDY, varargin);

STUDY = default_params(STUDY);

function STUDY = default_params(STUDY)
    if ~isfield(STUDY.etc, 'erpimparams'), STUDY.etc.erpimparams = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'erpimageopt'),  STUDY.etc.erpimparams.erpimageopt  = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'sorttype'   ),  STUDY.etc.erpimparams.sorttype  = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'sortwin'    ),  STUDY.etc.erpimparams.sortwin   = []; end;
    if ~isfield(STUDY.etc.erpimparams, 'sortfield'  ),  STUDY.etc.erpimparams.sortfield = 'latency'; end;
    if ~isfield(STUDY.etc.erpimparams, 'statistics' ),  STUDY.etc.erpimparams.statistics = 'param'; end;
    if ~isfield(STUDY.etc.erpimparams, 'groupstats' ),  STUDY.etc.erpimparams.groupstats = 'off'; end;
    if ~isfield(STUDY.etc.erpimparams, 'condstats'  ),  STUDY.etc.erpimparams.condstats  = 'off'; end;
    if ~isfield(STUDY.etc.erpimparams, 'threshold'  ),  STUDY.etc.erpimparams.threshold = NaN; end;
    if ~isfield(STUDY.etc.erpimparams, 'naccu') ,       STUDY.etc.erpimparams.naccu     = []; end;

