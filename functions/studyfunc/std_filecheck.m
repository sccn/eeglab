% std_filecheck() - Check if ERSP or SPEC file contain specific parameters.
%                  This file must contain a Matlab structure with a field 
%                  named 'parameter'. The content of this field will be
%                  compared to the 'params' input. If they are identical
%                  the output flag will indicate that recomputing this 
%                  file is not necessary. If they are different, the user
%                  is queried ('guion' option) to see if he wishes to use
%                  the new parameters and recompute the file (not done in
%                  this function) or if he wishes to use the parameters
%                  of the file on disk.
%
% Usage:  
% >> [ recompflag params ] = std_filecheck(filename, params, mode, ignorefields);
%
% Inputs:
%   filename     - [string] file containing a given measure (ERSP data for 
%                  instance). 
%   params       - [cell array or structure] cell array of parameters or
%                  structure. This is compared to the 'parameters' field 
%                  in the data file.
%   mode         - ['guion'|'usedisk'|'recompute'] 'guion' query the user 
%                  if the disk and input parameters are different. The
%                  outcome may be either 'usedisk' or 'recompute'. See
%                  recompflag output for more information.
%   ignorefields - [cell array] list fields to ignore
%
% Outputs:
%   recompflag - ['same'|'different'|'recompute'|'usedisk'|'cancel'] 'same' 
%                (resp. 'different') indicates that the parameters in the 
%                data file are identical (resp. different). 'recompute' 
%                indicate that the measure should be recomputed and the 
%                file has been erased. 'usedisk' indicate that the user
%                wishes (from the GUI) to use the version on disk.
%   params     - [structure] final parameter. This is usually identical
%                to the 'params' input unless the user choose to use
%                parameters from the file on disk. These are then
%                copied to this output structure.
%
% Authors: Arnaud Delorme, 2006-

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

function [ res, params2 ] = std_filecheck(filename, params2, guiflag, ignorefields);
    
    if nargin < 2
        help std_filecheck;
        return;
    end;
    if nargin < 3
        guiflag = 'guion';
    end;
    if nargin < 4
        ignorefields = {};
    end;
    
    if ~exist( filename ), res = guiflag; return; end;
        
    params1 = load('-mat', filename, 'parameters');
    params1 = finputcheck( params1.parameters, { 'tmp' 'real' [] NaN}, '', 'ignore'); % allow to tackle duplicate fields
    params1 = rmfield(params1, 'tmp');
    if iscell(params2), params2 = struct(params2{:}); end;
    
    % test if the fields are different
    % --------------------------------
    params1 = orderfields(params1);
    params2 = orderfields(params2);
    fields1 = fieldnames( params1 ); fields1 = setdiff( fields1, ignorefields);
    fields2 = fieldnames( params2 ); fields2 = setdiff( fields2, ignorefields);
    allfields = union(fields1, fields2);
    
    % make fields the same
    % --------------------
    res = 'same';
    if ~isequal( fields1, fields2 ),
        for ind = 1:length(allfields)
            if strcmpi(allfields{ind}, 'plotitc'), adsfads; end;
            if ~isfield( params1, allfields{ind})
                if isstr(getfield(params2, allfields{ind}))
                     params1 = setfield(params1, allfields{ind}, '');
                else params1 = setfield(params1, allfields{ind}, []);
                end;
                res = 'different';
            end;
            if ~isfield( params2, allfields{ind})
                if isstr(getfield(params1, allfields{ind}))
                     params2 = setfield(params2, allfields{ind}, '');
                else params2 = setfield(params2, allfields{ind}, []);
                end;
                res = 'different';
            end;
        end;
    end;
    
    % data type
    % ---------
    if ~isempty(strmatch('cycles', allfields)), strcom = 'ERSP';
    else                                        strcom = 'SPECTRAL';
    end;
    
    % compare fields
    % --------------
    txt = {};
    for ind = 1:length(allfields)
        if ~strcmp(allfields{ind},'baseline')
            val1 = getfield(params1, allfields{ind});
            val2 = getfield(params2, allfields{ind});
            val1str = fastif(isempty(val1), 'not set', vararg2str({val1(1:min(3, length(val1)))}));
            val2str = fastif(isempty(val2), 'not set', vararg2str({val2(1:min(3, length(val2)))}));

            tmptxt = sprintf('        ''%s'' is %s in the file (vs. %s)', allfields{ind}, val1str, val2str);
            if length(val1) ~= length(val2)
                res        = 'different';
                txt{end+1} = tmptxt;
            elseif ~isequal(val1, val2)
                if ~isnan(val1) & ~isnan(val2)
                    res        = 'different';
                    txt{end+1} = tmptxt;
                end
            end
        end
    end
            
    % build gui or return
    % -------------------
    if strcmpi(guiflag, 'usedisk'),
        if isempty(txt)
            params2 = params1;
            res     = 'usedisk';
            disp(['Using file on disk: ' filename ]);
            return; 
        else
            strvcat(txt{:})
            error([ 'Two ' strcom ' files had different parameters and the ' strcom ' function' 10 ...
                          'cannot handle that. We suggest that you delete all files. See command line details' ]);
        end;
    elseif strcmpi(guiflag, 'recompute'), 
        res     = 'recompute';
        disp(['Deleting and recomputing file: ' filename ]);
        return;
    elseif strcmpi(res, 'same') & ( strcmpi(guiflag, 'guion') | strcmpi(guiflag, 'same') )
        disp(['Using file on disk: ' filename ]);
        return;
    end;
    
    set_yes = [ 'set(findobj(''parent'', gcbf, ''tag'', ''ersp_no''), ''value'', 0);'];
    set_no  = [ 'set(findobj(''parent'', gcbf, ''tag'', ''ersp_yes''), ''value'', 0);' ];
    textgui1 = strvcat( [  upper(strcom)  ' info exists in file: ' filename ], ...
                'However, as detailed below, it does not fit with the requested values:');
    textgui2 = strvcat(txt{:});
    
    %textgui2 = ['wavelet cycles - [' num2str(params1.cycles(1)) ' ' num2str(params2.cycles(2)) ...
    %                '] instead of [' num2str(params1.cycles(1)) ' ' num2str(params2.cycles(2)) ...
    %            '], padratio - ' num2str(params1.padratio) ' instead of ' num2str(param2.padratio) ...
    %            ', and bootstrap significance - ' num2str(params1.alpha) ' instead of ' num2str(params2.alpha) ];
    
    uilist = { {'style' 'text' 'string' textgui1 } ...
               {'style' 'text' 'string' textgui2 } {} ...
               {'style' 'text' 'string' ['Would you like to recompute ' upper(strcom) ' and overwrite those values?' ]} ...
               {'style' 'checkbox' 'tag' 'ersp_yes' 'string' 'Yes, recompute' 'value' 1 'Callback' set_yes }  ...
               {'style' 'checkbox' 'tag' 'ersp_no' 'string' 'No, use data file parameters (and existing ERSP info)' 'value' 0 'Callback' set_no } };
    
    ersp_ans = inputgui('geometry', {[1] [1] [1] [1] [0.5 1] }, 'geomvert', [2 max(length(txt),1)*0.7 1 1 1 1], 'uilist', uilist, ...
                        'helpcom', '', 'title', ['Recalculate ' upper(strcom) ' parameters -- part of std_ersp()']); 
                        
    if isempty(ersp_ans), res = 'cancel'; return; end;
    if find(celltomat(ersp_ans))  == 2 % use existing ERSP info from this dataset
        disp(['Using file on disk: ' filename ]);
        params2 = params1;
        res              = 'usedisk';
    else % Over write data in dataset
        disp(['Deleting and recomputing file: ' filename ]);
        res              = 'recompute';
        %delete(filename);
    end
    
