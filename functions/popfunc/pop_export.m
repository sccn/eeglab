% pop_export() - export EEG dataset
%
% Usage:
%   >> com = pop_export(EEG);   % a window pops up
%   >> com = pop_export(EEG, filename, 'key', 'val', ... );
%
% Inputs:
%   EEG            - eeglab dataset
%   filename       - file name
%
% Optional inputs:
%   'ica'          - ['on'|'off'] export ICA activities (or ERP). Default 'off'.
%   'time'         - ['on'|'off'] include time values. Default 'on'.
%   'timeunit'     - [float] time unit rel. to seconds. Default: 1E-3 = msec.
%   'elec'         - ['on'|'off'] include electrodes names or component numbers. 
%                    Default 'on'.
%   'transpose'    - ['on'|'off'] 'off'-> electrode data = rows; 'on' -> electrode
%                    data = columns. Default 'off'.
%   'erp'          - ['on'|'off'] export ERP instead of raw data. Default 'off'.
%   'expr'         - [string] evaluate epxression on data. The expression must 
%                    contain a variable 'x' representing the 2-D or 3-D
%                    data. For example "x = 2*x" to multiply the data by 2.
%   'precision'    - [float] number of significant digits in output. Default 7.
%                    Default of 7 should allow to reach about 23 to 24 bits
%                    of precision and should be enough for EEG.
% 
% Outputs:
%   com            - The expresion that execute this function. i.e. 'pop_export(MyEEG, 'ExpEEG.mat')'
%
% Note: tabulation are used as a delimiter.
%
% Author: Arnaud Delorme, CNL / Salk Institute, May 13, 2003

% Copyright (C) May 13, 2003, Arnaud Delorme, Salk Institute, arno@salk.edu
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


function com = pop_export(EEG, filename, varargin); 

com = '';
if nargin < 1 
    help pop_export;
    return;
end

if nargin < 2
    commandload = [ '[filename, filepath] = uiputfile(''*'', ''Select a text file'');' ...
                    'if filename ~=0,' ...
                    '   set(findobj(''parent'', gcbf, ''tag'', ''tagedit''), ''string'', [ filepath filename ]);' ...
                    'end;' ...
                    'clear filename filepath tagtest;' ];
    
   uilist = { { 'style' 'text' 'string' 'Output file name' }, ...
              { 'style' 'edit' 'string' '' 'tag' 'tagedit' }, ...
              { 'style' 'pushbutton' 'string' 'Browse' 'callback' commandload },...
              { 'style' 'text' 'string' 'Export ICA activities instead of EEG data:' }, ...
              { 'style' 'checkbox' 'string' '' }, { },{ },  ...
              { 'style' 'text' 'string' 'Export ERP average instead of trials:' }, ...
              { 'style' 'checkbox' 'string' '' }, { }, { }, ...
              { 'style' 'text' 'string' 'Transpose matrix (elec -> cols):' }, ...
              { 'style' 'checkbox' 'string' '' 'value' 1}, { }, { }, ...
              { 'style' 'text' 'string' 'Export channel labels/component numbers:' }, ...
              { 'style' 'checkbox' 'string' '' 'value' 1 }, { }, { }, ...
              { 'style' 'text' 'string' 'Use comma separator (csv) instead of tab:' }, ...
              { 'style' 'checkbox' 'string' '' 'value' 0 }, { }, { }, ...
              { 'style' 'text' 'string' 'Export time values:' }, ...
              { 'style' 'checkbox' 'string' '' 'value' 1 }, ...
              { 'style' 'text' 'string' 'Unit (re. sec)' }, { 'style' 'edit' 'string' '1E-3' }, ...
              { 'style' 'text' 'string' 'Number of significant digits to output:' }, ...
              { 'style' 'edit' 'string' '4' }, { }, { }, ...
              { 'style' 'text' 'string' 'Apply an expression to the output (see ''expr'' help):'} , ...
              { 'style' 'edit' 'string' '' } { } { } };
   bcheck = [1.7 0.25 0.7 0.6];
   uigeom = { [1 2 0.5] bcheck bcheck bcheck bcheck bcheck bcheck [3 0.7 0.8 0.6] [3 2 0.05 0.05] };
   result = inputgui( uigeom, uilist, 'pophelp(''pop_export'');', 'Export data - pop_export()' );
   if isempty( result ), return; end
   
   % decode options
   % --------------
   if isempty(result{1}), error('File name required'); end
   filename = result{1};
   options = {}; 
   if result{2},  options = { options{:} 'ica' 'on' }; end
   if result{3},  options = { options{:} 'erp' 'on' }; end
   if result{4},  options = { options{:} 'transpose' 'on' }; end
   if ~result{5}, options = { options{:} 'elec' 'off' }; end
   if result{6},  options = { options{:} 'separator' ',' }; end
   if ~result{7}, options = { options{:} 'time' 'off' }; end
   if ~strcmpi(result{8}, '1E-3'), options = { options{:} 'timeunit' eval(result{8}) }; end
   if ~strcmpi(result{9}, '7'),    options = { options{:} 'precision' eval(result{9}) }; end
   if ~isempty(result{10}), options = { options{:} 'expr' result{10} }; end
else
    options = varargin;
end

% test inputs
% -----------
g = finputcheck(options, { ...
    'ica'       'string'    { 'on';'off' }     'off';
    'time'      'string'    { 'on';'off' }     'on';
    'timeunit'  'float'     [0 Inf]            1E-3;
    'elec'      'string'    { 'on';'off' }     'on';
    'transpose' 'string'    { 'on';'off' }     'off';
    'erp'       'string'    { 'on';'off' }     'off';
    'precision' 'integer'   [0 Inf]            7;
    'separator' 'string'    []                 char(9);
    'expr'      'string'    []                 '' }, 'pop_export');
if ischar(g), error(g); end

% select data
% ----------
if strcmpi(g.ica, 'on')
    eeglab_options; % changed from eeglaboptions 3/30/02 -sm
	if option_computeica  
		x = EEG.icaact;
	else
        x = EEG.icaweights*EEG.icasphere*reshape(EEG.data(EEG.icachansind,:,:), length(EEG.icachansind), EEG.trials*EEG.pnts);
        x = reshape(x, size(x,1), EEG.pnts, EEG.trials);
    end
else
    x = EEG.data;
end

% select erp
% ----------
if strcmpi(g.erp, 'on')
    x = mean(x, 3);
else 
    x = reshape(x, size(x,1), size(x,2)*size(x,3));
end

% write data
% ----------
if ~isempty(g.expr)
    eval([ g.expr ';' ]);
end

% add time axis
% -------------
if strcmpi(g.time, 'on')
    timeind = repmat( linspace(EEG.xmin, EEG.xmax, EEG.pnts)/g.timeunit, ...
                      [ 1 fastif(strcmpi(g.erp,'on'), 1, EEG.trials) ]);    
    xx = zeros(size(x,1)+1, size(x,2));
    xx(1,:)     = timeind;
    xx(2:end,:) = x;
    x = xx; clear xx;
end

% transpose and write to disk
% ---------------------------
fid = fopen(filename, 'w');
strprintf = [ '%.' num2str(g.precision) 'f' g.separator ];
if strcmpi(g.transpose, 'on')
    % columns
    % -------
    if strcmpi(g.elec, 'on')
        if strcmpi(g.time, 'on')
            fprintf(fid, 'Time%c', g.separator);
        end
        for index = 1:EEG.nbchan
            if ~isempty(EEG.chanlocs) && ~strcmpi(g.ica, 'on')
                fprintf(fid, '%s', EEG.chanlocs(index).labels);
            else fprintf(fid, '%d', index);
            end
            if index ~= EEG.nbchan
                fprintf(fid, '%c', g.separator);
            else
                fprintf(fid, '\n');
            end
        end
    end
    for index = 1:size(x,2)
        str = sprintf(strprintf, x(:,index));
        fprintf(fid, '%s\n', str(1:end-1));
    end
else 
    % writing electrodes
    % ------------------
    for index = 1:size(x,1)
        if strcmpi(g.time, 'on'), tmpind = index-1;
        else                      tmpind = index;
        end
        if strcmpi(g.elec, 'on') 
            if tmpind > 0
                if ~isempty(EEG.chanlocs) && ~strcmpi(g.ica, 'on')
                    fprintf(fid,'%s%c', EEG.chanlocs(tmpind).labels, g.separator );
                else fprintf(fid,'%d%c', tmpind, g.separator );
                end
            else 
                fprintf(fid, ' Time%c', g.separator);
            end
        end
        str = sprintf( [ '%.' num2str(g.precision) 'f' g.separator  ], x(index, :));
        fprintf(fid, '%s\n', str(1:end-1));
    end
end
fclose(fid);

com = sprintf('pop_export(EEG,%s);', vararg2str({ filename, options{:} })); 
return;
