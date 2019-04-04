% hlp_argstruct2linearcell() - Linearize configation output of arg_guipanel
%
% Usage:
%   >> cellval = hlp_argstruct2linearcells( cfg );
%
% Inputs:
%   cfg  - output configuration structure from arg_guipanel
%
% Output:
%   cellval - cell array of output values
%
% Author: Arnaud Delorme, SCCN & CERCO, CNRS, 2013-

% Copyright (C) 2013 Arnaud Delorme
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

function cellval = hlp_argstruct2linearcell(cfg);

    cellval = {};
    if isstruct(cfg)
        ff = fieldnames(cfg);
        for iField = 1:length(ff)
            if ~strcmpi(ff{iField}, 'arg_direct')
                if strcmpi(ff{iField}, 'arg_selection')
                    cellval = { cfg.arg_selection cellval{:} };
                else
                    val = cfg.(ff{iField});
                    if isstruct(val)
                        val = hlp_argstruct2linearcell(val);
                        cellval = { cellval{:} ff{iField} val{:} };
                    elseif iscell(val)
                        cellval = { cellval{:} ff{iField} vararg2str(val) };
                    else
                        cellval = { cellval{:} ff{iField} val };
                    end
                end
            end
        end
    else
        cellval = cfg;
    end
    
            
