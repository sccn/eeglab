% pop_rejchanspec() - reject artifacts channels in an EEG dataset using 
%                  channel spectrum. The average spectrum for all selected
%                  is computed and a threshold is applied.
%
% Usage:
%   >> pop_rejchanspec( INEEG ) % pop-up interactive window mode
%   >> [OUTEEG, indelec] = pop_rejchanspec( INEEG, 'key', 'val');
%
% Inputs:
%   INEEG       - input EEGLAB dataset
%
% Optional inputs:
%   'freqlims'  - [min max] frequency limits. May also be an array where
%                 each row defines a different set of limits. Default is 
%                 35 to the Niquist frequency of the data.
%   'stdthresh' - [max] positive threshold in terms of standard deviation.
%                 Default is 5.
%   'absthresh' - [max] positive threshold in terms of spectrum units
%                 (overrides the option above).
%   'averef'    - ['on'|'off'] 'on' computes average reference before
%                 applying threshold. Default is 'off'.
%   'plothist'  - ['on'|'off'] 'on' plot the histogram of values along 
%                 with the threshold.
%   'plotchans'  - ['on'|'off'] 'on' plot the channels scrollplot with
%                 selected channels for rejection in red. Allow selected
%                 channels rejection with the 'REJECT' button.
%   'elec'      - [integer array] only include specific channels. Default
%                 is to use all channels.
%   'specdata'  - [fload array] use this array containing the precomputed 
%                 spectrum instead of computing the spectrum. Default is
%                 empty.
%   'specfreqs' - [fload array] frequency array for precomputed spectrum
%                 above.
%   'verbose'   - ['on'|'off'] display information. Default is 'off'.
%
% Outputs:
%   OUTEEG    - output dataset with updated joint probability array
%   indelec   - indices of rejected electrodes
%   specdata  - data spectrum for the selected channels
%   specfreqs - frequency array for spectrum above
%
% Author: Arnaud Delorme, CERCO, UPS/CNRS, 2008-

% Copyright (C) 2008 Arnaud Delorme, CERCO, UPS/CNRS
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

function [EEG allrmchan specdata specfreqs com] = pop_rejchanspec(EEG, varargin)

if nargin < 1
    help pop_rejchanspec;
    return;
end
allrmchan = [];
specdata  = [];
specfreqs = [];
com       = '';
if nargin < 2
    uilist = { { 'style' 'text' 'string' 'Electrode (number(s); Ex: 2 4 5)' } ...
               { 'style' 'edit' 'string' ['1:' int2str(EEG.nbchan)] } ...
               { 'style' 'text' 'string' 'Frequency limits [min max]' } ...
               { 'style' 'edit' 'string' [ '35 ' int2str(floor(EEG.srate/2)) ] } ...
               { 'style' 'text' 'string' 'Standard dev. threshold limits [max]' } ...
               { 'style' 'edit' 'string' '5' } ...
               { 'style' 'text' 'string' 'OR absolute threshold limit [min max]' } ...
               { 'style' 'edit' 'string' '' } ...
               { 'style' 'text' 'string' 'Compute average reference first (check=on)' } ...
               { 'style' 'checkbox' 'string' '' 'value' 0 } { } ...
               { 'style' 'text' 'string' 'Plot histogram of power values (check=on)' } ...
               { 'style' 'checkbox' 'string' '' 'value' 0 } { } ...               
               { 'style' 'text' 'string' 'Plot channels scrollplot (check=on)' } ...
               { 'style' 'checkbox' 'string' '' 'value' 0 } { } ...
             };
          
           
    geom = { [2 1] [2 1] [2 1] [2 1] [2 0.3 0.7] [2 0.3 0.7] [2 0.3 0.7] };
    result = inputgui( 'uilist', uilist, 'geometry', geom, 'title', 'Reject channel using spectrum -- pop_rejchanspec()', ...
        'helpcom', 'pophelp(''pop_rejchan'')');
    if isempty(result), return; end
    
    options = { 'elec' eval( [ '[' result{1} ']' ] ) 'stdthresh' str2num(result{3}) 'freqlims' str2num(result{2}) };
    if ~isempty(result{4})
        options = { options{:} 'absthresh' str2num(result{4}) };
    end
    if result{5}, 
         options = { options{:} 'averef', 'on' }; 
    end
    if result{6}, 
         options = { options{:} 'plothist', 'on' }; 
    end
    % Begin: Added by Romain on 22 July 2010
    if result{7}, 
         options = { options{:} 'plotchans', 'on' }; 
    end
    % End: Added by Romain on 22 July 2010
    
else
    options = varargin;
end

% decode options
% --------------
opt = finputcheck( options, { 'averef'    'string'    { 'on';'off' }       'off';
                              'plothist'  'string'    { 'on';'off' }       'off';
                              'plotchans' 'string'    { 'on';'off' }       'off';
                              'verbose'   'string'    { 'on';'off' }       'off';
                              'elec'      'integer'   []                   [1:EEG.nbchan];
                              'freqlims'  'real'   []                      [35 EEG.srate/2];
                              'specdata'  'real'   []                      [];
                              'specfreqs' 'real'   []                      [];
                              'absthresh' 'real'   []                      [];
                              'stdthresh' 'real'   []                      5 }, 'pop_rejchanspec');
if ischar(opt), error(opt); end

% compute average reference if necessary
if strcmpi(opt.averef, 'on')
     NEWEEG = pop_reref(EEG, [], 'exclude', setdiff([1:EEG.nbchan], opt.elec));
else NEWEEG = EEG;
end
if isempty(opt.specdata)
    [tmpspecdata specfreqs] = pop_spectopo(NEWEEG, 1, [], 'EEG' , 'percent', 100, 'freqrange',[0 EEG.srate/2], 'plot', 'off');
    % add back 0 channels
    devStd = std(EEG.data(:,:), [], 2);
    if any(devStd == 0)
        goodchan  = find(devStd ~= 0);
        specdata  = zeros(length(opt.elec), size(tmpspecdata,2));
        specdata(goodchan,:) = tmpspecdata;
    else
        specdata = tmpspecdata;
    end
else
    specdata  = opt.specdata;
    specfreqs = opt.specfreqs;
end

if size(opt.stdthresh,1) == 1 && size(opt.freqlims,1) > 1
    opt.stdthresh = ones(length(opt.stdthresh), size(opt.freqlims,1))*opt.stdthresh;  
end

allrmchan = [];
for index = 1:size(opt.freqlims,1)
    % select frequencies, compute median and std then reject channels
    % ---------------------------------------------------------------
    [tmp fbeg] = min(abs(specfreqs - opt.freqlims(index,1)));
    [tmp fend] = min(abs(specfreqs - opt.freqlims(index,2)));
    selectedspec = mean(specdata(opt.elec, fbeg:fend), 2);
    if ~isempty(opt.absthresh)
        rmchan = find(selectedspec <= opt.absthresh(1) | selectedspec >= opt.absthresh(2));
    else
        m = median(selectedspec);
        s = std( selectedspec);
        nbTresh = size(opt.stdthresh);
        if length(opt.stdthresh) > 1
            rmchan = find(selectedspec <= m+s*opt.stdthresh(index,1) | selectedspec >= m+s*opt.stdthresh(index,2));
        else 
            rmchan = find(selectedspec > m+s*opt.stdthresh(index));
        end
    end
    
    % print out results
    % -----------------
    if isempty(rmchan)
         textout = sprintf('Range %2.1f-%2.1f Hz: no channel removed\n',  opt.freqlims(index,1), opt.freqlims(index,2));
    else textout = sprintf('Range %2.1f-%2.1f Hz: channels %s removed\n', opt.freqlims(index,1), opt.freqlims(index,2), int2str(opt.elec(rmchan')));
    end
    fprintf(textout);
    if strcmpi(opt.verbose, 'on')
        for inde = 1:length(opt.elec)
            if ismember(inde, rmchan)
                 fprintf('Elec %s power: %1.2f *\n', EEG.chanlocs(opt.elec(inde)).labels, selectedspec(inde));
            else fprintf('Elec %s power: %1.2f\n', EEG.chanlocs(opt.elec(inde)).labels  , selectedspec(inde));
            end
        end
    end
    allrmchan = [ allrmchan rmchan' ];    
    
    % plot histogram
    % --------------
    if strcmpi(opt.plothist, 'on')
        figure; hist(selectedspec);
        hold on; yl = ylim;
        if ~isempty(opt.absthresh)   
            plot([opt.absthresh(1) opt.absthresh(1)], yl, 'r');
            plot([opt.absthresh(2) opt.absthresh(2)], yl, 'r');
        else
            if length(opt.stdthresh) > 1
                threshold1 =  m+s*opt.stdthresh(index,1);
                threshold2 =  m+s*opt.stdthresh(index,2);
                plot([m m], yl, 'g');
                plot([threshold1 threshold1], yl, 'r');
                plot([threshold2 threshold2], yl, 'r');
            else
                threshold =  m+s*opt.stdthresh(index,1);
                plot([threshold threshold], yl, 'r');
            end
        end
        title(textout);
    end
    
end
allrmchan = unique_bc(allrmchan);

com = sprintf('EEG = pop_rejchan(EEG, %s);', vararg2str(options));
if strcmpi(opt.plotchans, 'on')   
    tmpcom = [ 'EEGTMP = pop_select(EEG, ''nochannel'', [' num2str(opt.elec(allrmchan)) ']);' ];
    tmpcom = [ tmpcom ...
            'LASTCOM = ' vararg2str(com) ';' ...
            '[ALLEEG EEG CURRENTSET tmpcom] = pop_newset(ALLEEG, EEGTMP, CURRENTSET);' ...
            '   if ~isempty(tmpcom),' ... 
            '     EEG = eegh(LASTCOM, EEG);' ...
            '     eegh(tmpcom);' ...
            '     eeglab(''redraw'');' ...
            '  end; clear EEGTMP tmpcom;' ];
 
    colors = cell(1,length(opt.elec)); colors(:) = { 'k' };
    colors(allrmchan) = { 'r' }; colors = colors(end:-1:1);
    fprintf('%d electrodes labeled for rejection\n', length(find(allrmchan)));
    tmpchanlocs = EEG.chanlocs;
    if ~isempty(EEG.chanlocs), tmplocs = EEG.chanlocs(opt.elec); tmpelec = { tmpchanlocs(opt.elec).labels }';
    else                       tmplocs = []; tmpelec = mattocell([1:EEG.nbchan]');
    end
    eegplot(EEG.data(opt.elec,:,:), 'srate', EEG.srate, 'title', 'Scroll component activities -- eegplot()', ...
        'limits', [EEG.xmin EEG.xmax]*1000, 'color', colors, 'eloc_file', tmplocs, 'command', tmpcom);
else
    EEG = pop_select(EEG, 'nochannel', opt.elec(allrmchan));
end

if nargin < 2
    allrmchan = sprintf('EEG = pop_rejchanspec(EEG, %s);', vararg2str(options));
end
