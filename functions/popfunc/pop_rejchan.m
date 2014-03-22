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
%   'measure'  - ['prob'|'kurt'|'spec'] compute probability 'prob', kurtosis 'kurt'
%                or spectrum 'spec' for each channel. Default is 'kurt'.
%   'norm'     - ['on'|'off'] normalize measure above (using trimmed 
%                normalization as described in the function jointprob()
%                and rejkurt(). Default is 'off'.
%   'precomp'  - [float array] use this array instead of computing the 'prob' 
%                or 'kurt' measures.
%   'freqrange' - [min max] frequency range for spectrum computation.
%                Default is 1 to sampling rate divided by 2. The average
%                of the log spectral power is computed over the frequency 
%                range of interest.
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
    cb_select = [ 'if get(gcbo, ''value'') == 3,' ...
                  '   set(findobj(gcbf, ''tag'', ''spec''), ''enable'', ''on'');' ...
                  'else,' ...
                  '   set(findobj(gcbf, ''tag'', ''spec''), ''enable'', ''off'');' ...
                  'end;' ];
    cb_norm = [ 'if get(gcbo, ''value''),' ...
                  '   set(findobj(gcbf, ''tag'', ''normlab''), ''string'', ''Z-score threshold [max] or [min max]'');' ...
                  'else,' ...
                  '   set(findobj(gcbf, ''tag'', ''normlab''), ''string'', ''Absolute threshold [max] or [min max]'');' ...
                  'end;' ];
    uilist = { { 'style' 'text' 'string' 'Electrode (number(s); Ex: 2 4 5)' } ...
               { 'style' 'edit' 'string' ['1:' int2str(EEG.nbchan)] } ...
               { 'style' 'text' 'string' 'Measure to use' } ...
               { 'style' 'popupmenu' 'string' 'Probability|Kurtosis|Spectrum' 'value' 2 'callback' cb_select } ...
               { 'style' 'text' 'string' 'Normalize measure (check=on)' } ...
               { 'style' 'checkbox' 'string' '' 'value' 1 'callback' cb_norm } { } ...
               { 'style' 'text' 'string' 'Z-score threshold [max] or [min max]' 'tag' 'normlab' } ...
               { 'style' 'edit' 'string' '5' } ...
               { 'style' 'text' 'string' 'Spectrum freq. range' 'enable' 'off' 'tag' 'spec' } ...
               { 'style' 'edit' 'string' '1 8'  'enable' 'off' 'tag' 'spec' } };
    geom = { [2 1.3] [2 1.3] [2 0.4 0.9] [2 1.3] [2 1.3] };
    result = inputgui( 'uilist', uilist, 'geometry', geom, 'title', 'Reject channel -- pop_rejchan()', ...
        'helpcom', 'pophelp(''pop_rejchan'')');
    if isempty(result), return; end;
    
    options = { 'elec' eval( [ '[' result{1} ']' ] ) 'threshold' str2num(result{4}) };
    if result{3}, 
         options = { options{:} 'norm', 'on' }; 
    else options = { options{:} 'norm', 'off' }; 
    end;
    
    if result{2} == 1,     options = { options{:} 'measure', 'prob' };
    elseif result{2} == 2, options = { options{:} 'measure', 'kurt' }; 
    else                   options = { options{:} 'measure', 'spec' }; 
    end;

else
    options = varargin;
end;

opt = finputcheck( options, { 'norm'      'string'    { 'on';'off' }       'off';
                              'measure'   'string'    { 'prob';'kurt';'spec' }    'kurt';
                              'precomp'   'real'      []                   [];
                              'freqrange' 'real'      []                   [1 EEG.srate/2];
                              'elec'      'integer'   []                   [1:EEG.nbchan];
                              'threshold' 'real'   []                      400 }, 'pop_rejchan');
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
elseif strcmpi(opt.measure, 'kurt')
    fprintf('Computing kurtosis for channels...\n');
    [ measure indelec ] = rejkurt( reshape(EEG.data(opt.elec,:,:), length(opt.elec), size(EEG.data,2)*size(EEG.data,3)), opt.threshold, opt.precomp, normval);
else
    fprintf('Computing spectrum for channels...\n');
    measure = pop_spectopo(EEG, 1, [], 'EEG' , 'freqrange', opt.freqrange, 'plot','off');

    % consider that data below 20 db has been filtered and remove it
    indFiltered = find(mean(measure) < -20);
    if ~isempty(indFiltered) && indFiltered(1) > 11, measure = measure(:,1:indFiltered(1)-10); disp('Removing spectrum data below -20dB (most likelly filtered out)'); end;
    meanSpec = mean(measure);
    stdSpec  = std( measure);
    
%     for indChan = 1:size(measure,1)
%         if any(measure(indChan,:) > meanSpec+stdSpec*opt.threshold), indelec(indChan) = 1; end;
%     end;
    if strcmpi(opt.norm, 'on')
        measure1  = max(bsxfun(@rdivide, bsxfun(@minus, measure, meanSpec), stdSpec),[],2);
        if length(opt.threshold) > 1
            measure2 = min(bsxfun(@rdivide, bsxfun(@minus, measure, meanSpec), stdSpec),[],2);
            indelec = measure2 < opt.threshold(1) | measure1 > opt.threshold(end);
            disp('Selecting minimum and maximum normalized power over the frequency range');
        else
            indelec = measure1 > opt.threshold(1);
            disp('Selecting maximum normalized power over the frequency range');
        end;
    else
        measure1 = max(measure,[],2);
        if length(opt.threshold) > 1
            measure2 = min(measure,[],2);
            indelec = measure2 < opt.threshold(1) | measure1 > opt.threshold(end);
            disp('Selecting minimum and maximum power over the frequency range');
        else
            indelec = measure > opt.threshold(1);
            disp('Selecting maximum power over the frequency range');
        end;
    end;
    measure = measure1;
end;
colors = cell(1,length(opt.elec)); colors(:) = { 'k' };
colors(find(indelec)) = { 'r' }; colors = colors(end:-1:1);
fprintf('%d electrodes labeled for rejection\n', length(find(indelec)));

% output variables
indelec = find(indelec)';
tmpchanlocs = EEG.chanlocs;
if ~isempty(EEG.chanlocs), tmplocs = EEG.chanlocs(opt.elec); tmpelec = { tmpchanlocs(opt.elec).labels }';
else                       tmplocs = []; tmpelec = mattocell([1:EEG.nbchan]');
end;
if exist('measure2', 'var')
     fprintf('#\tElec.\t[min]\t[max]\n');
     tmpelec(:,3) = mattocell(measure2);
     tmpelec(:,4) = mattocell(measure);
else fprintf('#\tElec.\tMeasure\n');
     tmpelec(:,3) = mattocell(measure);
end;
tmpelec(:,2) = tmpelec(:,1);
tmpelec(:,1) = mattocell([1:length(measure)]');
for index = 1:size(tmpelec,1)
    if exist('measure2', 'var')
         fprintf('%d\t%s\t%3.2f\t%3.2f', tmpelec{index,1}, tmpelec{index,2}, tmpelec{index,3}, tmpelec{index,4});
    else fprintf('%d\t%s\t%3.2f'       , tmpelec{index,1}, tmpelec{index,2}, tmpelec{index,3});
    end;
    if any(indelec == index), fprintf('\t*Bad*\n');
    else                      fprintf('\n');
    end;
end;
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
    eegplot(EEG.data(opt.elec,:,:), 'srate', EEG.srate, 'title', 'Scroll component activities -- eegplot()', ...
			 'limits', [EEG.xmin EEG.xmax]*1000, 'color', colors(end:-1:1), 'eloc_file', tmplocs, 'command', tmpcom);
else
    EEG = pop_select(EEG, 'nochannel', opt.elec(indelec));
end;

return;
