% std_plottf() - plot ERSP/ITC images a component
%              or channel cluster in a STUDY. Also allows plotting scalp
%              maps.
% Usage:
%          >> std_plottf( times, freqs, data, 'key', 'val', ...)
% Inputs:
%  times - [vector] latencies in ms of the data points.
%  freqs - [vector] frequencies in Hz of the data points.
%  data  -  [cell array] mean data for each subject group and/or data
%           condition. For example, to plot mean ERPs from a STUDY 
%           for epochs of 800 frames in two conditions from three groups 
%           of 12 subjects:
%
%           >> data = { [800x12] [800x12] [800x12];... % 3 groups, cond 1
%                       [800x12] [800x12] [800x12] };  % 3 groups, cond 2
%           >> std_plottf(erp_ms,data);
%
%           By default, parametric statistics are computed across subjects 
%           in the three groups. (group,condition) ERP averages are plotted. 
%           See below and >> help statcond 
%           for more information about the statistical computations.
%
% Optional display parameters:
%  'datatype'    - ['ersp'|'itc'] data type {default: 'ersp'}
%  'titles'      - [cell array of string] titles for each of the subplots. 
%                  { default: none}
%
% Statistics options:
%  'groupstats'  - ['on'|'off'] Compute (or not) statistics across groups.
%                  {default: 'off'}
%  'condstats'   - ['on'|'off'] Compute (or not) statistics across groups.
%                  {default: 'off'}

%  'threshold'   - [NaN|real<<1] Significance threshold. NaN -> plot the 
%                  p-values themselves on a different figure. When possible, 
%                  significance regions are indicated below the data.
%                  {default: NaN}
%  'maskdata'    - ['on'|'off'] when threshold is non-NaN and not both 
%                  condition and group statistics are computed, the user 
%                  has the option to mask the data for significance.
%                  {defualt: 'off'}
%
% Other plotting options:
%  'plotmode'    - ['normal'|'condensed'] statistics plotting mode:
%                  'condensed' -> plot statistics under the curves 
%                  (when possible); 'normal' -> plot them in separate 
%                  axes {default: 'normal'}
%  'freqscale'   - ['log'|'linear'|'auto'] frequency plotting scale. This
%                  will only change the ordinate not interpolate the data.
%                  If you change this option blindly, your frequency scale
%                  might be innacurate {default: 'auto'}
%
% ITC/ERSP image plotting options:
%  'tftopoopt'   - [cell array] tftopo() plotting options (ERSP and ITC)
%  'caxis'       - [min max] color axis (ERSP, ITC, scalp maps)
%
% Scalp map plotting options:
%  'chanlocs'    - [struct] channel location structure
%
% Author: Arnaud Delorme, CERCO, CNRS, 2006-
%
% See also: pop_erspparams(), pop_erpparams(), pop_specparams(), statcond()

% Copyright (C) 2006 Arnaud Delorme
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

function [pgroup, pcond, pinter] = std_plottf(timevals, freqs, data, varargin)

pgroup = [];
pcond  = [];
pinter = [];
if nargin < 2
    help std_plottf;
    return;
end

opt = finputcheck( varargin, { 'titles'         'cell'   []              cellfun(@num2str, cell(20,20), 'uniformoutput', false);
                               'caxis'          'real'   []              [];
                               'ersplim'        'real'   []              []; % same as above
                               'itclim'         'real'   []              []; % same as above
                               'ylim'           'real'   []              [];
                               'tftopoopt'      'cell'   []              {};
                               'threshold'      'real'   []              NaN;
                               'unitx'          'string' []              'ms'; % just for titles
                               'unitcolor'      'string' {}              'dB';
                               'chanlocs'       'struct' []              struct('labels', {});
                               'freqscale'      'string' { 'log','linear','auto' }  'auto'; % note that paramsersp in std_erspplot contains the information as well
                               'effect'         'string' { 'main','marginal' }   'marginal';
                               'averagemode'    'string' { 'rms','ave' }   'rms';
                               'events'         'cell'   []              {};
                               'groupstats'     'cell'   []              {};
                               'condstats'      'cell'   []              {};
                               'interstats'     'cell'   []              {};                               
                               'maskdata'       'string' { 'on','off' }   'off';
                               'plottopo'       'string' { 'on','off' }   'off';
                               'datatype'       'string' { 'ersp','itc' 'erpim' }    'ersp';
                               'plotmode'       'string' { 'normal','condensed' }  'normal' }, 'std_plottf');
if ischar(opt), error(opt); end
if all(all(cellfun('size', data, 3)==1))               opt.singlesubject = 'on'; end

% remove empty entries
datapresent = ~cellfun(@isempty, data);
for c = size(data,1):-1:1, if sum(datapresent(c,:)) == 0, data(c,:) = []; opt.titles(c,:) = []; if ~isempty(opt.groupstats), opt.groupstats(c) = []; end; end; end
for g = size(data,2):-1:1, if sum(datapresent(:,g)) == 0, data(:,g) = []; opt.titles(:,g) = []; if ~isempty(opt.condstats ), opt.condstats( g) = []; end; end; end

if ~isempty(opt.groupstats) && ~isempty(opt.condstats) && strcmpi(opt.maskdata, 'on')
    disp('Cannot use ''maskdata'' option with both condition stat. and group stat. on');
    disp('Disabling statistics');
    opt.groupstats = {}; opt.condstats = {}; opt.maskdata = 'off'; 
end
if ~isempty(opt.ersplim), opt.caxis = opt.ersplim; end
if ~isempty(opt.itclim), opt.caxis = opt.itclim; end
onecol  = { 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' };
manycol = { 'b' 'r' 'g' 'k' 'c' 'y' };

nc = size(data,1);
ng = size(data,2);
if nc >= ng, opt.transpose = 'on';
else         opt.transpose = 'off';
end

% test log frequencies
% --------------------
if length(freqs) > 2 && strcmpi(opt.freqscale, 'auto')
    midind  = floor(length(freqs)/2);
    if abs(freqs(midind)/freqs(end) - 1/2) < 0.1, opt.freqscale = 'linear';
    else                                          opt.freqscale = 'log';
    end
end

% condensed plot
% --------------
if strcmpi(opt.plotmode, 'condensed') 
    meanplot = zeros(size(data{1},1), size(data{1},2));
    count = 0;
    for c = 1:nc
        for g = 1:ng
            if ~isempty(data{c,g})
                meanplot = meanplot + mean(data{c,g},3);
                count = count+1;
            end
        end
    end
    meanplot = meanplot/count;
    options = { 'chanlocs', opt.chanlocs, 'electrodes', 'off', 'cbar', 'on', ...
            'cmode', 'separate', opt.tftopoopt{:} };       
    if strcmpi(opt.datatype, 'erpim'), options = { options{:} 'ylabel' 'Trials' }; end
    if strcmpi(opt.freqscale, 'log'), options = { options{:} 'logfreq', 'native' }; end
    tftopo( meanplot', timevals, freqs, 'title', opt.titles{1}, options{:}); 
    currentHangle = gca;
    if ~isempty( opt.caxis )
        caxis( currentHangle, opt.caxis )
    end
    cbar_standard(opt.datatype, ng, opt.unitcolor); 
    axes(currentHangle); 
    return; 
end

% plotting paramters
% ------------------
if ng > 1 && ~isempty(opt.groupstats), addc = 1; else addc = 0; end
if nc > 1 && ~isempty(opt.condstats  ), addr = 1; else addr = 0; end

% compute significance mask
% --------------------------
pinterplot = {};
if strcmpi(opt.effect, 'marginal') || ng == 1 || nc == 1
    if ~isnan(opt.threshold) && ( ~isempty(opt.groupstats) || ~isempty(opt.condstats) )    
        pcondplot  = opt.condstats;
        pgroupplot = opt.groupstats;
        maxplot = 1;
    else
        for ind = 1:length(opt.condstats),  pcondplot{ind}  = -log10(opt.condstats{ind}); end
        for ind = 1:length(opt.groupstats), pgroupplot{ind} = -log10(opt.groupstats{ind}); end
        maxplot = 3;
    end
elseif strcmpi(opt.effect, 'main') && ~isempty(opt.interstats)
    if ~isnan(opt.threshold) && ( ~isempty(opt.groupstats) || ~isempty(opt.condstats) )    
        pcondplot  = { opt.interstats{1} };
        pgroupplot = { opt.interstats{2} };
        pinterplot = opt.interstats{3};
        maxplot = 1;
    else
        if ~isempty(opt.interstats{1}), pcondplot  = { -log10(opt.interstats{1}) }; end
        if ~isempty(opt.interstats{2}), pgroupplot = { -log10(opt.interstats{2}) }; end
        if ~isempty(opt.interstats{3}), pinterplot = -log10(opt.interstats{3}); end
        maxplot = 3;
    end
end

% -------------------------------
% masking for significance of not
% -------------------------------
statmask = 0;
if strcmpi(opt.maskdata, 'on') && ~isnan(opt.threshold) && ...
        (~isempty(opt.condstats) || ~isempty(opt.condstats))
    addc = 0; addr = 0; statmask = 1;
end

% -------------------------
% plot time/frequency image
% -------------------------
options = { 'chanlocs', opt.chanlocs, 'electrodes', 'off', 'cbar', 'off', ...
            'cmode', 'separate', opt.tftopoopt{:} };
if strcmpi(opt.freqscale, 'log'), options = { options{:} 'logfreq', 'native' }; end
if strcmpi(opt.datatype, 'erpim'), options = { options{:} 'ylabel' 'Trials' }; end
        
% adjust figure size
% ------------------
fig = figure('color', 'w');
pos = get(fig, 'position');
set(fig, 'position', [ pos(1)+15 pos(2)+15 pos(3)/2.5*(nc+addr), pos(4)/2*(ng+addc) ]);
pos = get(fig, 'position');
if strcmpi(opt.transpose, 'off'), set(gcf, 'position', [ pos(1) pos(2) pos(4) pos(3)]);
else                              set(gcf, 'position', pos);
end

% options
% -------
options = { 'limits' [NaN NaN NaN NaN opt.caxis] 'verbose' 'off' 'mode' opt.averagemode options{:} };

for c = 1:nc
    for g = 1:ng
        %hdl(c,g) = mysubplot(nc+addr, ng+addc, g + (c-1)*(ng+addc), opt.transpose);
        hdl(c,g) = mysubplot(nc+addr, ng+addc, c, g, opt.transpose);
        if ~isempty(data{c,g})
            if strcmpi(opt.plottopo, 'off')
                tmpplot = mean(data{c,g},3);
            else
                tmpplot = data{c,g};
                tmpplot = permute(tmpplot, [3 1 2 4]);
            end
            if ~isreal(tmpplot(1)), tmpplot = abs(tmpplot); end % comes second for processing single trials
            if statmask, 
                if ~isempty(opt.condstats), tmpplot(find(pcondplot{g}(:) == 0)) = 0;
                else                        tmpplot(find(pgroupplot{c}(:) == 0)) = 0;
                end
            end
            if ~isempty(opt.events) && ~isempty(opt.events{c,g})
                 tmpevents = mean(opt.events{c,g},2);
            else tmpevents = [];
            end
            if strcmpi(opt.plottopo, 'on') && length(opt.chanlocs) > 1
                metaplottopo(tmpplot, 'chanlocs', opt.chanlocs, 'plotfunc', 'tftopo', 'squeeze', 'on', ...
                    'plotargs', { timevals, freqs, 'events', tmpevents, options{:} }, 'title', opt.titles{c,g});
            else
                tftopo( tmpplot, timevals, freqs, 'events', tmpevents, 'title', opt.titles{c,g}, options{:});
            end

            if c > 1
                ylabel(''); 
            end
        end
    
        % statistics accross groups
        % -------------------------
        if strcmpi(opt.effect, 'marginal') || (strcmpi(opt.effect, 'main') && c == 1)
            if g == ng && ng > 1 && ~isempty(opt.groupstats) && ~isinf(pgroupplot{c}(1)) && ~statmask
                if strcmpi(opt.effect, 'main') && nc>1, centerc = nc/2-0.5; else centerc = 0; end
                hdl(c,g+1) = mysubplot(nc+addr, ng+addc, c+centerc, ng + 1, opt.transpose);
                pgroupplot{c}(pgroupplot{c}<0) = 0;
                tmpOptions = { 'limits' [nan nan nan nan -maxplot maxplot] options{3:end} };
                if strcmpi(opt.plottopo, 'on') && length(opt.chanlocs) > 1
                    metaplottopo(permute(pgroupplot{c}, [3 1 2]), 'chanlocs', opt.chanlocs, 'plotfunc', 'tftopo', 'squeeze', 'on', ...
                        'plotargs', { timevals, freqs, tmpOptions{:} }, 'title', opt.titles{c,g+1});
                else
                    tftopo( pgroupplot{c}, timevals, freqs, 'title', opt.titles{c,g+1}, tmpOptions{:});
                end
            end
        end
    end
end

for g = 1:ng
    % statistics accross conditions
    % -----------------------------
    if strcmpi(opt.effect, 'marginal') || (strcmpi(opt.effect, 'main') && g == 1)
        if ~isempty(opt.condstats) && ~isinf(pcondplot{g}(1)) && ~statmask && nc > 1
            if strcmpi(opt.effect, 'main') && ng>1, centerg = ng/2-0.5; else centerg = 0; end
            hdl(nc+1,g) = mysubplot(nc+addr, ng+addc, nc+addr, g+centerg, opt.transpose);
            pcondplot{g}(pcondplot{g}<0) = 0;
            tmpOptions = { 'limits' [nan nan nan nan -maxplot maxplot] options{3:end} };
            if strcmpi(opt.plottopo, 'on') && length(opt.chanlocs) > 1
                metaplottopo(permute(pcondplot{g}, [3 1 2]), 'chanlocs', opt.chanlocs, 'plotfunc', 'tftopo', 'squeeze', 'on', ...
                    'plotargs', { timevals, freqs, tmpOptions{:} }, 'title', opt.titles{nc+1,g});
            else
                tftopo( pcondplot{g}, timevals, freqs, 'title', opt.titles{nc+1,g}, options{:});
            end
        end
    end
end

% statistics accross group and conditions
% ---------------------------------------
if ~isempty(opt.groupstats) && ~isempty(opt.condstats) && ng > 1 && nc > 1 && ~isempty(pinterplot)
    hdl(nc+1,ng+1) = mysubplot(nc+addr, ng+addc, nc+addr, ng+1, opt.transpose);
    pinterplot(pinterplot<0) = 0;
    tftopo( pinterplot,  timevals, freqs, 'title', opt.titles{nc+1,ng+1}, options{:});
    caxis([-maxplot maxplot]);
    ylabel('');
end    

% color bars
% ----------
axes(hdl(nc,ng)); 
cbar_standard(opt.datatype, ng, opt.unitcolor); 
if isnan(opt.threshold) && (nc ~= size(hdl,1) || ng ~= size(hdl,2))
    ind = find(ishandle(hdl(end:-1:1)));
    axes(hdl(end-ind(1)+1));
    cbar_signif(ng, maxplot);
end

% mysubplot2 (allow to transpose if necessary)
% -------------------------------------------
function hdl = mysubplot(nr,nc,r,c,subplottype)

    cmargin = 0.2/nc;
    rmargin = 0.2/nr;
    if strcmpi(subplottype, 'transpose') || strcmpi(subplottype, 'on'),   hdl = subplot('position',[(r-1)/nr+rmargin (nc-c)/nc+cmargin 1/nr-2*rmargin 1/nc-2*cmargin]);
    elseif strcmpi(subplottype, 'normal') || strcmpi(subplottype, 'off'), hdl = subplot('position',[(c-1)/nc+cmargin (nr-r)/nr+rmargin 1/nc-2*cmargin 1/nr-2*rmargin]);
    elseif strcmpi(subplottype, 'noplot'), hdl = gca;
    else error('Unknown subplot type');
    end
    
% % mysubplot (allow to transpose if necessary)
% % -------------------------------------------
% function hdl = mysubplot(nr,nc,ind,transp);
% 
%     r = ceil(ind/nc);
%     c = ind -(r-1)*nc;
%     if strcmpi(transp, 'on'), hdl = subplot(nc,nr,(c-1)*nr+r);
%     else                      hdl = subplot(nr,nc,(r-1)*nc+c);
%     end

% colorbar for ERSP and scalp plot
% --------------------------------
function cbar_standard(datatype, ng, unitcolor);
    pos = get(gca, 'position');
    tmpc = caxis;
    fact = fastif(ng == 1, 40, 20);
    tmp = axes('position', [ pos(1)+pos(3)+max(pos(3)/fact,0.006) pos(2) max(pos(3)/fact,0.01) pos(4) ]);  
    set(gca, 'unit', 'normalized');
    if strcmpi(datatype, 'itc')
         cbar(tmp, 0, tmpc, 10); ylim([0.5 1]);
         title('ITC','fontsize',10,'fontweight','normal','interpreter','none');
    elseif strcmpi(datatype, 'erpim')
        cbar(tmp, 0, tmpc, 5);
    else
        cbar(tmp, 0, tmpc, 5);
        title(unitcolor);
    end
    

% colorbar for significance
% -------------------------
function cbar_signif(ng, maxplot);
    % Retrieving Defaults
    icadefs;
    
    pos = get(gca, 'position');
    tmpc = caxis;
    fact = fastif(ng == 1, 40, 20);
    tmp = axes('position', [ pos(1)+pos(3)+max(pos(3)/fact,0.006) pos(2) max(pos(3)/fact,0.01) pos(4) ]);  
    map = colormap(DEFAULT_COLORMAP);
    n = size(map,1);
    cols = [ceil(n/2):n]';
    image([0 1],linspace(0,maxplot,length(cols)),[cols cols]);
    %cbar(tmp, 0, tmpc, 5);
    tick = linspace(0, maxplot, maxplot+1);
    set(gca, 'ytickmode', 'manual', 'YAxisLocation', 'right', 'xtick', [], ...
        'ytick', tick, 'yticklabel', round(10.^-tick*1000)/1000);
    xlabel('');
    colormap(DEFAULT_COLORMAP);

% rapid filtering for ERP
% -----------------------
function tmpdata2 = myfilt(tmpdata, lowpass, highpass, factor, filtertype)

    tmpdata2 = reshape(tmpdata, size(tmpdata,1), size(tmpdata,2)*size(tmpdata,3)*size(tmpdata,4));
    tmpdata2 = eegfiltfft(tmpdata2',lowpass, highpass, factor, filtertype)';
    tmpdata2 = reshape(tmpdata2, size(tmpdata,1), size(tmpdata,2), size(tmpdata,3), size(tmpdata,4));
