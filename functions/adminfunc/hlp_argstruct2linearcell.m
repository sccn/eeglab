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
                    end;
                end;
            end;
        end
    else
        cellval = cfg;
    end;
    
            
