% pop_rejchan() - reject artifacts channels in an EEG dataset using joint 
%                  probability of the recorded electrode.
%
% Usage:
%   >> pop_rejchan( INEEG ) % pop-up interative window mode
%   >> [OUTEEG, locthresh, globthresh, nrej] = ...
%		= pop_rejchan( INEEG, 'key', 'val');
%
% Inputs:
%   INEEG      - input dataset
%
% Optional inputs:
%   'elec'     - [n1 n2 ...] electrode number(s) to take into 
%                consideration for rejection
%   'threshold' - [max] absolute thresold or activity probability 
%                 limit(s) (in std. dev.) if norm is 'on'.
%   'measure'  - ['prob'|'kurt'] compute probability 'prob' or kurtosis 'kurt'
%                for each channel. Default is 'kurt'.
%   'norm'     - ['on'|'off'] normalize measure above (using trimmed 
%                normalization as described in the function jointprob()
%                and rejkurt(). Default is 'off'.
%   'precomp'  - [float array] use this array instead of computing the 'prob' 
%                or 'kurt' measures.
%
% Outputs:
%   OUTEEG    - output dataset with updated joint probability array
%   indelec   - indices of rejected electrodes
%   measure   - measure value for each electrode
%
% Author: Arnaud Delorme, CERCO, UPS/CNRS, 2008-
%
% See also: jointprob(), rejkurt()

% Copyright (C) 2008 Arnaud Delorme, CERCO, UPS/CNRS
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

function [EEG, indelec, measure, com] = pop_rejchan( EEG, varargin);

com = '';
indelec = [];
measure = [];
if nargin < 1
   help pop_rejchan;
   return;
end;  

if nargin < 2

	% which set to save
	% -----------------
    uilist = { { 'style' 'text' 'string' 'Electrode (number(s); Ex: 2 4 5)' } ...
               { 'style' 'edit' 'string' ['1:' int2str(EEG.nbchan)] } ...
               { 'style' 'text' 'string' 'Measure to use' } ...
               { 'style' 'popupmenu' 'string' 'Probability|Kurtosis' 'value' 2 } ...
               { 'style' 'text' 'string' 'Normalize measure (check=on)' } ...
               { 'style' 'checkbox' 'string' '' 'value' 1 } { } ...
               { 'style' 'text' 'string' 'Threshold limits [max]' } ...
               { 'style' 'edit' 'string' '5' } };
    geom = { [2 1.3] [2 1.3] [2 0.4 0.9] [2 1.3] };
    result = inputgui( 'uilist', uilist, 'geometry', geom, 'title', 'Reject channel -- pop_rejchan()', ...
        'helpcom', 'pophelp(''pop_rejchan'')');
    if isempty(result), return; end;
    
    options = { 'elec' eval( [ '[' result{1} ']' ] ) 'threshold' str2num(result{4}) };
    if result{3}, 
         options = { options{:} 'norm', 'on' }; 
    else options = { options{:} 'norm', 'off' }; 
    end;
    
    if result{2} == 1, options = { options{:} 'measure', 'prob' };
    else               options = { options{:} 'measure', 'kurt' }; 
    end;

else
    options = varargin;
end;

opt = finputcheck( options, { 'norm'    'string'    { 'on';'off' }       'off';
                              'measure' 'string'    { 'prob';'kurt' }    'kurt';
                              'precomp' 'real'      []                   [];
                              'elec'    'integer'   []                   [1:EEG.nbchan];
                              'threshold' 'real'   []                    400 }, 'pop_rejchan');
if isstr(opt), error(opt); end;

% compute the joint probability
% -----------------------------
if strcmpi(opt.norm, 'on')
    normval = 2;
else
    normval = 0;
end;
if strcmpi(opt.measure, 'prob')
    fprintf('Computing probability for channels...\n');
    [ measure indelec ] = jointprob( reshape(EEG.data(opt.elec,:,:), length(opt.elec), size(EEG.data,2)*size(EEG.data,3)), opt.threshold, opt.precomp, normval);
else
    fprintf('Computing kurtosis for channels...\n');
    [ measure indelec ] = rejkurt( reshape(EEG.data(opt.elec,:,:), length(opt.elec), size(EEG.data,2)*size(EEG.data,3)), opt.threshold, opt.precomp, normval);
end;
colors = cell(1,length(opt.elec)); colors(:) = { 'k' };
colors(find(indelec)) = { 'r' }; colors = colors(end:-1:1);
fprintf('%d electrodes labeled for rejection\n', length(find(indelec)));

indelec = find(indelec)';
if isempty(indelec), return; end;

com = sprintf('EEG = pop_rejchan(EEG, %s);', vararg2str(options));
if nargin < 2
    tmpcom = [ 'EEGTMP = pop_select(EEG, ''nochannel'', [' num2str(opt.elec(indelec)) ']);' ];
    tmpcom = [ tmpcom ...
            'LASTCOM = ' vararg2str(com) ';' ...
            '[ALLEEG EEG CURRENTSET tmpcom] = pop_newset(ALLEEG, EEGTMP, CURRENTSET);' ...
            '   if ~isempty(tmpcom),' ... 
            '     EEG = eegh(LASTCOM, EEG);' ...
            '     eegh(tmpcom);' ...
            '     eeglab(''redraw'');' ...
            '  end; clear EEGTMP tmpcom;' ];
    tmpchanlocs = EEG.chanlocs;
    if ~isempty(EEG.chanlocs), tmplocs = EEG.chanlocs(opt.elec); tmpelec = { tmpchanlocs(opt.elec).labels }';
    else                       tmplocs = []; tmpelec = mattocell([1:EEG.nbchan]');
    end;
    eegplot(EEG.data(opt.elec,:,:), 'srate', EEG.srate, 'title', 'Scroll component activities -- eegplot()', ...
			 'limits', [EEG.xmin EEG.xmax]*1000, 'color', colors, 'eloc_file', tmplocs, 'command', tmpcom);
    
    tmpelec(:,3) = mattocell(measure);
    tmpelec(:,2) = tmpelec(:,1);
    tmpelec(:,1) = mattocell([1:length(measure)]')
else
    EEG = pop_select(EEG, 'nochannel', opt.elec(indelec));
end;

return;
