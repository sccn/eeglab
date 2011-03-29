% pop_multifit() - fit multiple component dipoles using DIPFIT 
%
% Usage:
%         >> EEG = pop_multifit(EEG); % pop-up graphical interface
%         >> EEG = pop_multifit(EEG, comps, 'key', 'val', ...);
%
% Inputs:
%  EEG      - input EEGLAB dataset.
%  comps    - indices component to fit. Empty is all components.
%
% Optional inputs:
%  'dipoles'   - [1|2] use either 1 dipole or 2 dipoles contrain in
%                symmetry. Default is 1.
%  'dipplot'   - ['on'|'off'] plot dipoles. Default is 'off'.
%  'plotopt'   - [cell array] dipplot() 'key', 'val' options. Default is
%                'normlen', 'on', 'image', 'fullmri'
%  'rmout'     - ['on'|'off'] remove dipoles outside the head. Artifactual
%                component often localize outside the head. Default is 'off'.
%  'threshold' - [float] rejection threshold during component scan.
%                Default is 40 (residual variance above 40%).
%
% Outputs:
%  EEG      - output dataset with updated "EEG.dipfit" field
%
% Note: residual variance is set to NaN if DIPFIT does not converge
%  
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, Oct. 2003

% Copyright (C) 9/2003 Arnaud Delorme, SCCN/INC/UCSD, arno@salk.edu
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

function [EEG, com] = pop_multifit(EEG, comps, varargin);
    
    if nargin < 1
        help pop_multifit;
        return;
    end;
    
    com = [];
    ncomps = size(EEG.icaweights,1);
    if ncomps == 0, error('you must run ICA first'); end;

    if nargin<2
        cb_chans = 'tmplocs = EEG.chanlocs; set(findobj(gcbf, ''tag'', ''chans''), ''string'', int2str(pop_chansel({tmplocs.labels}))); clear tmplocs;';
        
        uilist = { { 'style' 'text' 'string' 'Component indices' } ...
                   { 'style' 'edit' 'string' [ '1:' int2str(ncomps) ] } ... 
                    { 'style' 'text' 'string' 'Rejection threshold RV (%)' } ...
                   { 'style' 'edit' 'string' '100' } ... 
                   { 'style' 'text' 'string' 'Remove dipoles outside the head' } ...
                   { 'style' 'checkbox' 'string' '' 'value' 0 } {} ...
                   { 'style' 'text' 'string' 'Fit bilateral dipoles (check)' } ...
                   { 'style' 'checkbox' 'string' '' 'value' 0 } {} ...
                   { 'style' 'text' 'string' 'Plot resulting dipoles (check)' } ...
                   { 'style' 'checkbox' 'string' '' 'value' 0 } {} ...
                   { 'style' 'text' 'string' 'dipplot() plotting options' } ...
                   { 'style' 'edit' 'string' '''normlen'' ''on''' } ...
                   { 'style' 'pushbutton' 'string' 'Help' 'callback' 'pophelp(''dipplot'')' } }; 

        results = inputgui( { [1.91 2.8] [1.91 2.8] [3.1 0.8 1.6] [3.1 0.8 1.6] [3.1 0.8 1.6] [2.12 2.2 0.8]}, ...
                            uilist, 'pophelp(''pop_multifit'')', ...
                            'Fit multiple ICA components -- pop_multifit()');
        if length(results) == 0 return; end;
        comps        = eval( [ '[' results{1} ']' ] );
        
        % selecting model
        % ---------------
        options = {};
         if ~isempty(results{2})
            options      = { options{:} 'threshold' eval( results{2} ) };
        end;
        if results{3}, options = { options{:} 'rmout' 'on' }; end;
        if results{4}, options = { options{:} 'dipoles' 2 }; end;
        if results{5}, options = { options{:} 'dipplot' 'on' }; end;
        options = { options{:} 'plotopt' eval( [ '{ ' results{6} ' }' ]) };
    else 
        options = varargin;
    end;
    
    % checking parameters
    % -------------------
    if isempty(comps), comps = [1:size(EEG.icaweights,1)]; end;
    g = finputcheck(options, { 'settings'  { 'cell' 'struct' }    []        {};                % deprecated
                               'dipoles'   'integer'  [1 2]      1;
                               'threshold' 'float'    [0 100]   40;
                               'dipplot'   'string'   { 'on' 'off' } 'off';
                               'rmout'     'string'   { 'on' 'off' } 'off';
                               'plotopt'   'cell'     {}        {'normlen' 'on' }});
    
    if isstr(g), error(g); end;    
    EEG     = eeg_checkset(EEG, 'chanlocs_homogeneous');
    
    % dipfit settings
    % ---------------
    if isstruct(g.settings)
        EEG.dipfit = g.settings;
    elseif ~isempty(g.settings)
        EEG = pop_dipfit_settings( EEG, g.settings{:}); % will probably not work but who knows
    end;
            
    % Scanning dipole locations
    % -------------------------
    dipfitdefs;
    skipscan = 0;
    try 
        alls = cellfun('size', { EEG.dipfit.model.posxyz }, 2);
        if length(alls) == ncomps
            if all(alls == 3)
                skipscan = 1;
            end;
        end;
    catch, end;
    if skipscan
        disp('Skipping scanning since all dipoles have non-null starting positions.');
    else
        disp('Scanning dipolar grid to find acceptable starting positions...');
        xg  = linspace(-floor(meanradius), floor(meanradius),11);
        yg  = linspace(-floor(meanradius), floor(meanradius),11);
        zg  = linspace(0                 , floor(meanradius), 6);
        EEG = pop_dipfit_gridsearch( EEG, [1:ncomps], ...
                                eval(xgridstr), eval(ygridstr), eval(zgridstr), 100);
        disp('Scanning terminated. Refining dipole locations...');
    end;
   
    % set symmetry constraint
    % ----------------------
    if strcmpi(EEG.dipfit.coordformat,'MNI')
        defaultconstraint = 'x';
    else
        defaultconstraint = 'y';
    end;
    
    % Searching dipole localization
    % -----------------------------
    disp('Searching dipoles locations...');
    chansel =  EEG.dipfit.chansel;
    %elc     = getelecpos(EEG.chanlocs, EEG.dipfit);
    plotcomps = [];
    for i = comps(:)'
        if i <= length(EEG.dipfit.model) & ~isempty(EEG.dipfit.model(i).posxyz)
            if g.dipoles == 2,
                % try to find a good origin for automatic dipole localization
                EEG.dipfit.model(i).active = [1 2];
                EEG.dipfit.model(i).select = [1 2];
                if isempty(EEG.dipfit.model(i).posxyz)
                    EEG.dipfit.model(i).posxyz = zeros(1,3);
                    EEG.dipfit.model(i).momxyz = zeros(2,3);
                else
                    EEG.dipfit.model(i).posxyz(2,:) = EEG.dipfit.model(i).posxyz;
                    if strcmpi(EEG.dipfit.coordformat, 'MNI')
                         EEG.dipfit.model(i).posxyz(:,1) = [-40;40];
                    else EEG.dipfit.model(i).posxyz(:,2) = [-40;40];
                    end;
                    EEG.dipfit.model(i).momxyz(2,:) = EEG.dipfit.model(i).momxyz;
                end;
            else 
                EEG.dipfit.model(i).active = [1];
                EEG.dipfit.model(i).select = [1];
            end;
            warning backtrace off;
            try,
                if g.dipoles == 2,
                    EEG = dipfit_nonlinear(EEG, 'component', i, 'symmetry', defaultconstraint);
                else
                    EEG = dipfit_nonlinear(EEG, 'component', i, 'symmetry', []);
                end;
            catch, EEG.dipfit.model(i).rv = NaN; disp('Maximum number of iterations reached. Fitting failed');
            end;
            warning backtrace on;
            plotcomps = [ plotcomps i ];
        end;
    end;
    
    % set RV to 1 for dipole with higher than 40% residual variance
    % -------------------------------------------------------------
    EEG.dipfit.model  = dipfit_reject(EEG.dipfit.model, g.threshold/100);

    % removing dipoles outside the head
    % ---------------------------------
    if strcmpi(g.rmout, 'on') & strcmpi(EEG.dipfit.coordformat, 'spherical')
        rmdip = [];
        for index = plotcomps
            if ~isempty(EEG.dipfit.model(index).posxyz)
                if any(sqrt(sum(EEG.dipfit.model(index).posxyz.^2,2)) > 85)
                    rmdip = [ rmdip index];
                    EEG.dipfit.model(index).posxyz = [];
                    EEG.dipfit.model(index).momxyz = [];
                    EEG.dipfit.model(index).rv     = 1;
                end;
            end;
        end;
        plotcomps = setdiff(plotcomps, rmdip);
        if length(rmdip) > 0
            fprintf('%d out of cortex dipoles removed (usually artifacts)\n', length(rmdip));
        end;
    end;
    
    % plotting dipoles
    % ----------------
    if strcmpi(g.dipplot, 'on')
        pop_dipplot(EEG, 'DIPFIT', plotcomps, g.plotopt{:});
    end;
    
    com = sprintf('%s = pop_multifit(%s, %s);', inputname(1), inputname(1), vararg2str({ comps options{:}}));
    return;
    
% get electrode positions from eeglag
% -----------------------------------
function elc = getelecpos(chanlocs, dipfitstruct);
    try,
        elc = [ [chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]' ];
    catch
        disp('No 3-D carthesian coordinates; re-computing them from 2-D polar coordinates');
        EEG.chanlocs = convertlocs(EEG.chanlocs, 'topo2all');
        elc = [ [chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]' ];
    end;
    % constrain electrode to sphere
    % -----------------------------
    disp('Constraining electrodes to sphere');
    elc = elc - repmat( dipfitstruct.vol.o, [size(elc,1) 1]); % recenter
    % (note the step above is not needed since the origin should always be 0)
    elc = elc ./ repmat( sqrt(sum(elc.*elc,2)), [1 3]); % normalize
    elc = elc * max(dipfitstruct.vol.r);         % head size

    %for index= 1:size(elc,1)
    %    elc(index,:) = max(dipfitstruct.vol.r) * elc(index,:) /norm(elc(index,:));
    %end;


    
