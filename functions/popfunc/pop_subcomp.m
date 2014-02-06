% pop_subcomp() - remove specified components from an EEG dataset.
%                 and subtract their activities from the data. Else,
%                 remove components already marked for rejection.
% Usage:
%   >> OUTEEG = pop_subcomp( INEEG ); % pop-up window mode
%   >> OUTEEG = pop_subcomp( INEEG, components, confirm);
%
% Pop-up window interface:
%   "Component(s) to remove ..." - [edit box] Array of components to 
%                remove from the data. Sets the 'components' parameter 
%                in the command line call (see below).
%   "Component(s) to retain ..." - [edit box] Array of components to
%                to retain in the data. Sets the 'components' parameter in
%                the command line call. Then, comp_to_remove = ...
%                    setdiff([1:size(EEG.icaweights,1)], comp_to_keep)
%                Overwrites "Component(s) to remove" (above).
% Command line inputs:
%   INEEG      - Input EEG dataset.
%   components - Array of components to remove from the data. If empty, 
%                 remove components previously marked for rejection (e.g., 
%                 EEG.reject.gcompreject).
%   confirm    - [0|1] Display the difference between original and processed
%                dataset. 1 = Ask for confirmation. 0 = Do not ask. {Default: 0}
% Outputs:
%   OUTEEG     - output dataset.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: compvar()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 
% 02-15-02 propagate ica weight matrix -ad sm jorn 

function [EEG, com] = pop_subcomp( EEG, components, plotag )

com='';
if nargin < 1
   help pop_subcomp;
   return;
end;
if nargin < 3
	plotag = 0;
end;	
if nargin < 2
	% popup window parameters
	% -----------------------
	if ~isempty(EEG.reject.gcompreject)
        components = find(EEG.reject.gcompreject == 1);
        components = components(:)';
        promptstr    = { ['Component(s) to remove from the data ([] = marked comps.)'] };
        %promptstr    = { ['Components to subtract from data' 10 '(default: pre-labeled components to reject):'] };
    else
        components = [];
        promptstr    = { ['Component(s) to remove from data:'] };
    end;
    uilist    = { { 'style' 'text' 'string' ['Component(s) to remove from data:'] } ...
                  { 'style' 'edit' 'string' int2str(components) } ...
                  { 'style' 'text' 'string' 'Component(s) to retain (overwrites "Component(s) to remove")' } ...
                  { 'style' 'edit' 'string' '' } ...
                  };
    geom = { [2 0.7] [2 0.7] };
	result       = inputgui( 'uilist', uilist, 'geometry', geom, 'helpcom', 'pophelp(''pop_subcomp'')', ...
                                     'title', 'Remove components from data -- pop_subcomp()');
	if length(result) == 0 return; end;
	components   = eval( [ '[' result{1} ']' ] );
    if ~isempty(result{2}), 
        components   = eval( [ '[' result{2} ']' ] );
        components  = setdiff_bc([1:size(EEG.icaweights,1)], components);
    end;
end;
 
if isempty(components)
	if ~isempty(EEG.reject.gcompreject)
      		components = find(EEG.reject.gcompreject == 1);
   	else
        	fprintf('Warning: no components specified, no rejection performed\n');
         	return;
   	end;
else
    if (max(components) > size(EEG.icaweights,1)) || min(components) < 1
        error('Component index out of range');
    end;
end;

fprintf('Computing projection ....\n');
component_keep = setdiff_bc(1:size(EEG.icaweights,1), components);
compproj = EEG.icawinv(:, component_keep)*eeg_getdatact(EEG, 'component', component_keep, 'reshape', '2d');
compproj = reshape(compproj, size(compproj,1), EEG.pnts, EEG.trials);

%fprintf( 'The ICA projection accounts for %2.2f percent of the data\n', 100*varegg);
	
if nargin < 2 | plotag ~= 0

    ButtonName = 'continue';
    while ~strcmpi(ButtonName, 'Cancel') & ~strcmpi(ButtonName, 'Accept')
        ButtonName=questdlg2( [ 'Please confirm. Are you sure you want to remove these components?' ], ...
                             'Confirmation', 'Cancel', 'Plot ERPs', 'Plot single trials', 'Accept', 'Accept');
        if strcmpi(ButtonName, 'Plot ERPs')
            if EEG.trials > 1
                tracing  = [ squeeze(mean(EEG.data(EEG.icachansind,:,:),3)) squeeze(mean(compproj,3))];
                figure;   
                plotdata(tracing, EEG.pnts, [EEG.xmin*1000 EEG.xmax*1000 0 0], ...
                    'Trial ERPs (red) with and (blue) without these components');
            else
                warndlg2('Cannot plot ERPs for continuous data');
            end;
        elseif strcmpi(ButtonName, 'Plot single trials')  
        	eegplot( EEG.data(EEG.icachansind,:,:), 'srate', EEG.srate, 'title', 'Black = channel before rejection; red = after rejection -- eegplot()', ...
            	 'limits', [EEG.xmin EEG.xmax]*1000, 'data2', compproj); 
        end;
    end;    
    switch ButtonName,
        case 'Cancel', 
        	disp('Operation cancelled');
        	return; 
        case 'Accept',
       		disp('Components removed');
    end % switch
end;
EEG.data(EEG.icachansind,:,:) = compproj;
EEG.setname = [ EEG.setname ' pruned with ICA'];
EEG.icaact  = [];
goodinds    = setdiff_bc(1:size(EEG.icaweights,1), components);
EEG.icawinv     = EEG.icawinv(:,goodinds);
EEG.icaweights  = EEG.icaweights(goodinds,:);
EEG.specicaact  = [];
EEG.specdata    = [];
EEG.reject      = [];

try,
    EEG.dipfit.model = EEG.dipfit.model(goodinds);
catch, end;

com = sprintf('%s = pop_subcomp( %s, [%s], %d);', inputname(1), inputname(1), ...
   int2str(components), plotag);
return;
