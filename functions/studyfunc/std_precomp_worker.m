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
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function feature = std_precomp_worker(filename, varargin)

if nargin < 2
    help std_precomp_worker;
    return;
end

g = struct(varargin{2:end});

% load dataset
% ------------
[STUDY ALLEEG] = pop_loadstudy('filename', filename);
if ~isfield(g, 'design'), g.design = STUDY.currentdesign; end
if isfield(g, 'erp') || isfield(g, 'spec') || isfield(g, 'erpim') || isfield(g, 'scalp')
    error('This function is currently designed to compute time-frequency decompositions only'); 
end
if ~isfield(g, 'ersp') && ~isfield(g, 'itc')
    error('You must compute either ERSP or ITC when using the EC2 cluster'); 
end

% run std_precomp (THIS IS THE PART WE WANT TO PARALELIZE)
% ---------------
% for index = 1:length(STUDY.design(g.design).cell)
%     [STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, varargin{:}, 'cell', index);
% end
std_precomp(STUDY, ALLEEG, varargin{:});

filebase = STUDY.design(g.design).cell(g.cell).filebase;
if isfield(g, 'channel')
    fileERSP = [ filebase '.datersp' ];
    fileITC  = [ filebase '.datitc' ];
    if exist(fileERSP), feature.ersp = load('-mat', fileERSP); end
    if exist(fileITC ), feature.itc  = load('-mat', fileITC ); end
else
    fileERSP = [ filebase '.icaersp' ];
    fileITC  = [ filebase '.icaitc' ];
    if exist(fileERSP), feature.ersp = load('-mat', fileERSP); end
    if exist(fileITC ), feature.itc  = load('-mat', fileITC ); end
end
