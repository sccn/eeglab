% std_pacplot() - Commandline function to plot cluster PAC 
%                (phase-amplitude coupling).
%
% Usage:    
%   >> [STUDY] = std_pacplot(STUDY, ALLEEG, key1, val1, key2, val2);  
%   >> [STUDY itcdata itctimes itcfreqs pgroup pcond pinter] = ...
%                std_pacplot(STUDY, ALLEEG ...);
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the datasets in the STUDY. 
%                Note: ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%
% Additional help:
% Inputs and output of this function are strictly identical to the std_erspplot(). 
% See the help message of this function for more information.
%
% See also: std_erspplot(), pop_clustedit(), pop_preclust()
%
% Authors: Arnaud Delorme, CERCO, July, 2009-

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

function [STUDY, allpac, alltimes, allfreqs, pgroup, pcond, pinter] = std_pacplot(STUDY, ALLEEG, varargin)

if nargin < 2
    help std_pacplot;
    return;
end

[STUDY allpac alltimes allfreqs pgroup pcond pinter ] = std_erspplot(STUDY, ALLEEG, 'datatype', 'pac', varargin{:});
