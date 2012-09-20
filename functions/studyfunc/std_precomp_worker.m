% std_precomp_worker() - allow dispatching ERSP to be computed in parallel 
%                        on a given cluster. 
% Usage:
%   >> feature = std_precomp_worker(filename, varargin);
%
% Inputs:
%   filename - STUDY file name
%
% Optional inputs:
%   Optional inputs are the same as for the std_precomp function. Note that 
%   this function can currently only compute ERSP and ITC. The argument
%   'cell' must be defined so a given node will only compute the measure on
%   one cell (one cell per node).
%                   output trials {default: whole measure range}
% Output:
%  feature  - data structure containing ERSP and/or ITC
%
% Author: Arnaud Delorme, SCCN, UCSD, 2012-

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

function feature = std_precomp_worker(filename, varargin)

if nargin < 2
    help std_precomp_worker;
    return;
end;

g = struct(varargin{2:end});

% load dataset
% ------------
[STUDY ALLEEG] = pop_loadstudy('filename', filename);
if ~isfield(g, 'design'), g.design = STUDY.currentdesign; end;
if isfield(g, 'erp') || isfield(g, 'spec') || isfield(g, 'erpim') || isfield(g, 'scalp')
    error('This function is currently designed to compute time-frequency decompositions only'); 
end;
if ~isfield(g, 'ersp') && ~isfield(g, 'itc')
    error('You must compute either ERSP or ITC when using the EC2 cluster'); 
end;

% run std_precomp (THIS IS THE PART WE WANT TO PARALELIZE)
% ---------------
% for index = 1:length(STUDY.design(g.design).cell)
%     [STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, varargin{:}, 'cell', index);
% end;
std_precomp(STUDY, ALLEEG, varargin{:});

filebase = STUDY.design(g.design).cell(g.cell).filebase;
if isfield(g, 'channel')
    fileERSP = [ filebase '.datersp' ];
    fileITC  = [ filebase '.datitc' ];
    if exist(fileERSP), feature.ersp = load('-mat', fileERSP); end;
    if exist(fileITC ), feature.itc  = load('-mat', fileITC ); end;
else
    fileERSP = [ filebase '.icaersp' ];
    fileITC  = [ filebase '.icaitc' ];
    if exist(fileERSP), feature.ersp = load('-mat', fileERSP); end;
    if exist(fileITC ), feature.itc  = load('-mat', fileITC ); end;
end;
