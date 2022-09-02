% pop_subcomp() - remove specified components from an EEG dataset.
%                 and subtract their activities from the data. Else,
%                 remove components already marked for rejection. When used
%                 with the options 'keepcomp', the function will retain [1]
%                 or reject[0] the components provided as input.
% Usage:
%   >> OUTEEG = pop_subcomp( INEEG ); % pop-up window mode
%   >> OUTEEG = pop_subcomp( INEEG, components, plotag);
%   >> OUTEEG = pop_subcomp( INEEG, components, plotag, keepcomp);
%
% Pop-up window interface:
%   "Component(s) to remove ..." - [edit box] Array of components to 
%                remove from the data. Sets the 'components' parameter 
%                in the command line call (see below).
%   "Component(s) to retain ..." - [edit box] Array of components to
%                to retain in the data. Sets the 'components' parameter in
%                the command line call. Then, comp_to_remove = ...
%                    setdiff([1:size(EEG.icaweights,1)], comp_to_keep). See
%                    option 'keepcomp' for command line call.
%                Overwrites "Component(s) to remove" (above).
% Command line inputs:
%   INEEG      - Input EEG dataset.
%   components - Array of components to remove from the data. If empty, 
%                 remove components previously marked for rejection (e.g., 
%                 EEG.reject.gcompreject).
%   plotag     - [0|1] Display the difference between original and processed
%                dataset. 1 = Ask for confirmation. 0 = Do not ask. {Default: 0}
%   keepcomp   - [0|1] If [1] will retain the components provided by the
%                input variable 'components', [0] will reject them. Option
%                intended to be used only with the components provided in 
%                'components', and not with components marked for rejection
%                {Default: 0}
% Outputs:
%   OUTEEG     - output dataset.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: compvar()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad 
% 02-15-02 propagate ica weight matrix -ad sm jorn 

function [EEG, com] = pop_subcomp( EEG, components, plotag, keepflag)

com='';
if nargin < 1
   help pop_subcomp;
   return;
end
if nargin < 3
	plotag = 0;
end
if nargin < 4
    keepflag = 0;
end

if nargin < 2
    components = [];
	% popup window parameters
	% -----------------------
    res = 'Manual rej.';
    if any(cellfun(@(x)any(x.gcompreject), { EEG.reject }))
        if length(EEG) == 1
            compStr = sprintf('%d,', find(EEG.reject.gcompreject == 1));
            msg = sprintf('Components [%s] flagged for rejection.', compStr(1:end-1));
            res = questdlg2(strvcat(msg, 'Do you want to remove these components?', 'Note: we recommend removing components in STUDY instead'), 'Remove components from data', 'Cancel', 'Manual rej.', 'Yes', 'Cancel');
        else
            msg = 'Components flagged for rejection detected in some datasets.';
            res = questdlg2(strvcat(msg, 'Do you want to remove these components?', 'Note: we recommend removing components in STUDY instead'), 'Remove components from data', 'Cancel', 'Yes', 'Cancel');
        end
        if strcmpi(res, 'Cancel')
            return;
        end
        if strcmpi(res, 'Yes')
            components = '';
        end
    else
        if length(EEG) > 1
            warndlg2(strvcat('You have multiple datasets selected and no components', ...
                'flagged for rejection. Flag components first.'));
            return;
        end
    end
    
    if strcmpi(res, 'Manual rej.')
        if ~isempty(EEG.reject.gcompreject)
            components = find(EEG.reject.gcompreject == 1);
            components = components(:)';
            %promptstr    = { ['Components to subtract from data' 10 '(default: pre-labeled components to reject):'] };
        else
            components = [];
        end
        uilist    = { { 'style' 'text' 'string' 'Note: for group level analysis, remove components in STUDY' } ...
                      { 'style' 'text' 'string' 'List of component(s) to remove from data' } ...
                      { 'style' 'edit' 'string' int2str(components) } ...
                      { 'style' 'text' 'string' 'Or list of component(s) to retain' } ...
                      { 'style' 'edit' 'string' '' } ...
                      };
        geom = { 1 [2 0.7] [2 0.7] };
        result       = inputgui( 'uilist', uilist, 'geometry', geom, 'helpcom', 'pophelp(''pop_subcomp'')', ...
                                         'title', 'Remove components from data -- pop_subcomp()');
        if length(result) == 0 return; end
        components   = eval( [ '[' result{1} ']' ] );
        if ~isempty(result{2})
            componentsOld = components;
            components   = eval( [ '[' result{2} ']' ] );
            if isequal(components, componentsOld)
                components = [];
            end
            keepflag = 1; %components  = setdiff_bc([1:size(EEG.icaweights,1)], components); 
        end
    end
end
 
% process multiple datasets
% -------------------------
if length(EEG) > 1
    if nargin < 2
        [ EEG, com ] = eeg_eval( 'pop_subcomp', EEG, 'params', { components, plotag, keepflag }, 'warning', 'on' );
    else
        [ EEG, com ] = eeg_eval( 'pop_subcomp', EEG, 'params', { components, plotag, keepflag } );
    end
    if isempty( components )
        com = [ com ' % [] or '' means removing components flagged for rejection' ];
    end
    return;
end

componentsOri = components;
if isempty(components)
    if ~isempty(EEG.reject.gcompreject)
        components = find(EEG.reject.gcompreject == 1);
    else
        fprintf('Warning: no components specified, no rejection performed\n');
        return;
    end
else
    if keepflag == 1; components  = setdiff_bc([1:size(EEG.icaweights,1)], components); end
    if (max(components) > size(EEG.icaweights,1)) || min(components) < 1
        error('Component index out of range');
    end
end

fprintf('Computing projection and removing %d components ....\n', length(components));
component_keep = setdiff_bc(1:size(EEG.icaweights,1), components);
compproj = EEG.icawinv(:, component_keep)*eeg_getdatact(EEG, 'component', component_keep, 'reshape', '2d');
compproj = reshape(compproj, size(compproj,1), EEG.pnts, EEG.trials);

%fprintf( 'The ICA projection accounts for %2.2f percent of the data\n', 100*varegg);
	
if nargin < 2 || plotag ~= 0

    ButtonName = 'continue';
    while ~strcmpi(ButtonName, 'Cancel') && ~strcmpi(ButtonName, 'Accept')
        ButtonName=questdlg2( [ 'Please confirm your choice. Are you sure you want to remove the selected components from the data?' ], ...
                             'Confirmation', 'Cancel', 'Plot ERPs', 'Plot single trials', 'Accept', 'Accept');
        if strcmpi(ButtonName, 'Plot ERPs')
            if EEG.trials > 1
                tracing  = [ squeeze(mean(compproj,3)) squeeze(mean(EEG.data(EEG.icachansind,:,:),3))];
                figure;   
                plotdata(tracing, EEG.pnts, [EEG.xmin*1000 EEG.xmax*1000 0 0], ...
                    'Trial ERPs (blue) with and (red) without these components');
            else
                warndlg2('Cannot plot ERPs for continuous data');
            end
        elseif strcmpi(ButtonName, 'Plot single trials')  
        	eegplot( EEG.data(EEG.icachansind,:,:), 'srate', EEG.srate, 'title', 'Black = channel before rejection; red = after rejection -- eegplot()', ...
            	 'limits', [EEG.xmin EEG.xmax]*1000, 'data2', compproj); 
        end
    end
    switch ButtonName
        case 'Cancel'
        	disp('Operation cancelled');
        	return; 
        case 'Accept'
       		disp('Components removed');
    end % switch
end
EEG.data(EEG.icachansind,:,:) = compproj;
EEG.setname = [ EEG.setname ' pruned with ICA'];
EEG.icaact  = [];
goodinds    = setdiff_bc(1:size(EEG.icaweights,1), components);
EEG.icawinv     = EEG.icawinv(:,goodinds);
EEG.icaweights  = EEG.icaweights(goodinds,:);
EEG.specicaact  = [];
EEG.specdata    = [];
EEG.reject      = [];

if isfield(EEG.etc, 'ic_classification')
    if isfield(EEG.etc.ic_classification, 'ICLabel') 
        if isfield(EEG.etc.ic_classification.ICLabel, 'classifications')
            if ~isempty(EEG.etc.ic_classification.ICLabel.classifications)
                EEG.etc.ic_classification.ICLabel.classifications = EEG.etc.ic_classification.ICLabel.classifications(goodinds,:);
            end
        end
    end
end

try
    EEG.dipfit.model = EEG.dipfit.model(goodinds);
catch, end

com = sprintf('EEG = pop_subcomp( EEG, [%s], %d);', int2str(componentsOri(:)'), plotag);
if isempty( components )
    com = [ com ' % [] means removing components flagged for rejection' ];
end
return;
