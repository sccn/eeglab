% pop_multifit() - fit multiple component dipoles using DIPFIT 
%
% Usage:
%         >> EEG = pop_multifit(EEG); % pop-up graphical interface
%         >> EEG = pop_multifit(EEG, comps, 'key', 'val', ...);
%
% Inputs:
%  EEG      - input dataset.
%  comps    - component to fit
%
% Optional inputs:
%  'settings'  - [cell array] dipfit settings (arguments to the 
%                pop_dipfit_settings() function). Default is none.
%  'dipoles'   - [1|2] use either 1 dipole or 2 dipoles contrain in
%                symetry. Default is 1.
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

%123456789012345678901234567890123456789012345678901234567890123456789012

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

% $Log: not supported by cvs2svn $
% Revision 1.20  2004/05/20 18:06:46  arno
% header
%
% Revision 1.19  2004/03/04 19:28:37  arno
% email
%
% Revision 1.18  2003/12/16 00:22:31  arno
% header typo
%
% Revision 1.17  2003/12/05 23:14:23  arno
% msd
% msg
%
% Revision 1.16  2003/12/05 23:03:44  arno
% removing out of cortex dipoles
%
% Revision 1.15  2003/12/04 03:01:30  arno
% adding threshold parameter
%
% Revision 1.14  2003/11/16 18:15:55  scott
% Changing button '...' to 'Help' -sm
%
% Revision 1.13  2003/11/06 16:26:13  arno
% debug if # of ICA comps different from # of channells
%
% Revision 1.12  2003/11/06 16:01:00  arno
% header
%
% Revision 1.11  2003/11/05 16:22:52  arno
% genoushomogenous -> homogeneoushomogenous -> homogeneous
%
% Revision 1.10  2003/10/31 17:33:37  arno
% adding help buttons
%
% Revision 1.9  2003/10/31 17:26:58  arno
% adding electrode selection
% (I mean channel section)
%
% Revision 1.8  2003/10/31 01:31:08  scott
% text
%
% Revision 1.7  2003/10/30 01:18:23  arno
% scanning
%
% Revision 1.6  2003/10/29 23:50:49  arno
% scanning
%
% Revision 1.5  2003/10/29 16:28:13  arno
% allowing ploting and plotting options
%
% Revision 1.4  2003/10/29 15:13:43  arno
% adapting to new dipfit_settings
%
% Revision 1.3  2003/10/11 19:33:46  arno
% allowing fitting constrained in simetry
%
% Revision 1.2  2003/10/11 02:01:36  arno
% [Afixing history
%
% Revision 1.1  2003/10/11 01:37:05  arno
% Initial revision

function [EEG, com] = pop_multifit(EEG, comps, varargin);
    
    if nargin < 1
        help pop_multifit;
        return;
    end;
    
    com = [];
    ncomps = size(EEG.icaweights,1);
    if ncomps == 0, error('you must run ICA first'); end;

    if nargin<2
        cb_chans = 'set(findobj(gcbf, ''tag'', ''chans''), ''string'', int2str(pop_chansel(EEG.chanlocs)));';
        
        uilist = { { 'style' 'text' 'string' 'Component indices' } ...
                   { 'style' 'edit' 'string' [ '1:' int2str(ncomps) ] } ... 
                   { 'style' 'text' 'string' 'pop_dipfit_settings() options' } ...
                   { 'style' 'edit' 'string' '' } ... 
                   { 'style' 'pushbutton' 'string' 'Help' 'callback' 'pophelp(''pop_dipfit_settings'')' } ... 
                   { 'style' 'text' 'string' 'Omit the following channels' } ...
                   { 'style' 'edit' 'string' '' 'tag' 'chans' } ... 
                   { 'style' 'pushbutton' 'string' '...' 'callback' cb_chans } ... 
                   { 'style' 'text' 'string' 'Rejection threshold RV (%)' } ...
                   { 'style' 'edit' 'string' '40' } ... 
                   { 'style' 'text' 'string' 'Remove dipoles outside the head' } ...
                   { 'style' 'checkbox' 'string' '' 'value' 0 } {} ...
                   { 'style' 'text' 'string' 'Fit bilateral dipoles (check)' } ...
                   { 'style' 'checkbox' 'string' '' 'value' 0 } {} ...
                   { 'style' 'text' 'string' 'Plot resulting dipoles (check)' } ...
                   { 'style' 'checkbox' 'string' '' 'value' 0 } {} ...
                   { 'style' 'text' 'string' 'dipplot() plotting options' } ...
                   { 'style' 'edit' 'string' '''normlen'' ''on'' ''image'' ''fullmri''' } ...
                   { 'style' 'pushbutton' 'string' 'Help' 'callback' 'pophelp(''dipplot'')' } }; 

        results = inputgui( { [1.91 2.8] [2.12 2.2 0.8]  [2.12 2.2 0.8] [1.91 2.8] [3.1 0.4 2] [3.1 0.4 2] [3.1 0.4 2] [2.12 2.2 0.8]}, ...
                            uilist, 'pophelp(''pop_multifit'')', ...
                            'Fit multiple ICA components -- pop_multifit()');
        if length(results) == 0 return; end;
        comps        = eval( [ '[' results{1} ']' ] );
        if ~isempty(results{3})
             options      = { 'settings' eval( [ '{' results{2} ', ''electrodes'', setdiff(1:EEG.nbchan,[' results{3} ']) }' ] ) };
        else options      = { 'settings' eval( [ '{' results{2} '}' ] ) };
        end;
        if ~isempty(results{4})
            options      = { options{:} 'threshold' eval( results{4} ) };
        end;
        if results{5}, options = { options{:} 'rmout' 'on' }; end;
        if results{6}, options = { options{:} 'dipoles' 2 }; end;
        if results{7}, options = { options{:} 'dipplot' 'on' }; end;
        options = { options{:} 'plotopt' eval( [ '{ ' results{8} ' }' ]) };
    else 
        options = varargin;
    end;
    
    % checking parameters
    % -------------------
    g = finputcheck(options, { 'settings'  'cell'     []        {}; 
                               'dipoles'   'integer'  [1 2]      1;
                               'threshold' 'float'    [0 100]   40;
                               'dipplot'   'string'   { 'on' 'off' } 'off';
                               'rmout'     'string'   { 'on' 'off' } 'on';
                               'plotopt'   'cell'     {}        {'normlen' 'on' 'image' 'fullmri'}});
    
    if isstr(g), error(g); end;    
    EEG     = eeg_checkset(EEG, 'chanlocs_homogeneous');
    
    % removing eye channels
    % ---------------------
    %if length(EEG.chanlocs) > 69 & length(g.dipfit) < 3
    %[tmp tmp2 indreye] = intersect( 'reye', lower({ EEG.chanlocs.labels }));
    %[tmp tmp2 indleye] = intersect( 'leye', lower({ EEG.chanlocs.labels }));
    %elecsel = setdiff([1:71], [indleye indreye]);
    %disp('Removing eye channels');
    %elecopt = 0;
    %for index = 1:2:length(g.settings)
    %    if strcmpi(g.settings{index}, 'electrodes')
    %        g.settings{index+1} = intersect( g.settings{index+1}, elecsel);
    %        electopt = 1;
    %    end;
    %end;
    %if electopt
    %   EEG     = pop_dipfit_settings( EEG, g.settings{:});
    %else
    %    EEG     = pop_dipfit_settings( EEG, g.settings{:}, 'electrodes', elecsel);
    %end;
    if isempty(g.settings)
        EEG     = pop_dipfit_settings( EEG, 'electrodes', [1:EEG.nbchan]);        
    else
        EEG     = pop_dipfit_settings( EEG, g.settings{:} );        
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
        EEG = pop_dipfit_batch( EEG, [1:ncomps], ...
                                eval(xgridstr), eval(ygridstr), eval(zgridstr), 100);
        disp('Scanning terminated. Refining dipole locations...');
    end;
    EEG.dipfit.model  = dipfit_reject(EEG.dipfit.model, g.threshold/100);

    
    % Searching dipole localization
    % -----------------------------
    disp('Searching dipoles locations...');
    chansel =  EEG.dipfit.chansel;
    elc     = getelecpos(EEG.chanlocs, EEG.dipfit);
    plotcomps = [];
    for i = comps(:)'
        if i <= length(EEG.dipfit.model) & ~isempty(EEG.dipfit.model(i).posxyz)
            if g.dipoles == 2,
                EEG.dipfit.model(i).active = [1 2];
                EEG.dipfit.model(i).select = [1 2];
                EEG.dipfit.model(i).posxyz = zeros(2,3);
                EEG.dipfit.model(i).momxyz = zeros(2,3);
            else 
                EEG.dipfit.model(i).active = [1];
                EEG.dipfit.model(i).select = [1];
            end;
            try,
                if g.dipoles == 2,
                    EEG.dipfit.model(i) = dipfit_manual(EEG.dipfit.model(i), EEG.icawinv(chansel,i), ...
                                                        elc(chansel,:), EEG.dipfit.vol, 'method', 'position', 'constraint', defaultconstraint);
                else
                    EEG.dipfit.model(i) = dipfit_manual(EEG.dipfit.model(i), EEG.icawinv(chansel,i), ...
                                                        elc(chansel,:), EEG.dipfit.vol, 'method', 'position' );
                end;
            catch, EEG.dipfit.model(i).rv =NaN; disp('Maximum number of iterations reached. Fitting failed');
            end;
            plotcomps = [ plotcomps i ];
        end;
    end;
    
    % removing dipoles outside the head
    % ---------------------------------
    if strcmpi(g.rmout, 'on')
        rmdip = [];
        for index = plotcomps
            if ~isempty(EEG.dipfit.model(index).posxyz)
                if any(sqrt(sum(EEG.dipfit.model(index).posxyz.^2,2)) > EEG.dipfit.vol.r(1))
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
        elc = cell2mat({chanlocs.X; chanlocs.Y; chanlocs.Z}');
    catch
        disp('No 3-D carthesian coordinates; re-computing them from 2-D polar coordinates');
        EEG.chanlocs = convertlocs(EEG.chanlocs, 'topo2all');
        elc = cell2mat({chanlocs.X; chanlocs.Y; chanlocs.Z}');
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


    
