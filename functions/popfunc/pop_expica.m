% pop_expica() - export ICA weights or inverse matrix
%
% Usage:
%   >> pop_expica( EEG, whichica);             % a window pops up
%   >> pop_expica( EEG, whichica, filename );
%
% Inputs:
%   EEG         - EEGLAB dataset
%   whichica    - ['weights'|'inv'] export ica 'weights' or ica inverse
%                 matrix ('inv'). Note: for 'weights', the function 
%                 export the product of the sphere and weights matrix.
%   filename    - text file name
% 
% Author: Arnaud Delorme, CNL / Salk Institute, Mai 14, 2003
%
% See also: pop_export()

% Copyright (C) Mai 14, 2003, Arnaud Delorme, Salk Institute, arno@salk.edu
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

function com = pop_expica(EEG, whichica, filename); 
    
com = '';
if nargin < 1 
    help pop_expica;
    return;
end

if nargin < 2
    whichica = 'weights';
end
switch lower(whichica)
 case {'weights' 'inv'}, ;
 otherwise error('Unrecognized option for ''whichica'' parameter');
end

if nargin < 3
	% ask user
	[filename, filepath] = uiputfile('*.*', [ 'File name for ' ...
                        fastif(strcmpi(whichica, 'inv'), 'inverse', 'weight') ' matrix -- pop_expica()']); 
    drawnow;
	if filename == 0 return; end
	filename = [filepath filename];
end

% save datas
% ----------
if strcmpi(whichica, 'inv')
    tmpmat = double(EEG.icawinv);
else
    tmpmat = double(EEG.icaweights*EEG.icasphere);
end
save(filename, '-ascii', 'tmpmat');

com = sprintf('pop_expica(EEG, ''%s'', ''%s'');', whichica, filename); 

return;
