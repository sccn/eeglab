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
%  'channels'    - [cell array] channel names (for titles) {default: all}
%  'condnames'   - [cell array] names of conditions (for titles} 
%                  {default: none}
%  'groupnames'  - [cell array] names of subject groups (for titles)
%                  {default: none}
%  'subject'     - [string] plot subject name (for title)
%  'compinds'    - [integer] plot component index (for title)
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
%  'freqscale'   - ['log'|'linear'|'auto'] frequency plotting scale.
%                  {default: 'auto'}
%  'ylim'        - [min max] ordinate limits for ERP and spectrum plots
%                  {default: all available data}
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

% $Log: not supported by cvs2svn $
% Revision 1.18  2009/11/18 02:17:35  dev
% fix bug 782
%
% Revision 1.17  2009/10/07 05:07:19  arno
% Fix missing conditions/groups
%
% Revision 1.16  2009/08/29 04:24:56  arno
% new statistics
%
% Revision 1.15  2009/05/31 02:22:10  arno
% Adding FDR and bootstrap to all STUDY functions
%
% Revision 1.14  2007/09/11 10:50:30  arno
% fix numerous small display problems and crash
%
% Revision 1.13  2007/08/14 02:27:01  allen
% Same
%
% Revision 1.12  2007/08/14 02:26:16  allen
% fix data masking problem
%
% Revision 1.11  2007/08/13 21:06:33  nima
% _
%
% Revision 1.10  2007/08/12 23:53:32  arno
% new call to std_stat
%
% Revision 1.9  2007/08/10 19:47:40  nima
% _
%
% Revision 1.8  2007/08/07 01:06:09  allen
% Nima: Changed cbar for ERSP and P values, also removed redundant labels from y axis (Frequencies)
%
% Revision 1.7  2007/07/30 23:03:40  arno
% do not plot stats automatically
%
% Revision 1.6  2007/05/11 03:11:36  toby
% bug when plotting 'threshold' with more than one dataset fixed
%
% Revision 1.5  2007/04/06 01:18:16  arno
% frequency scaling
%
% Revision 1.4  2007/04/05 21:38:23  arno
% now plot in linear scale if necessaru
%
% Revision 1.3  2007/04/05 21:35:47  arno
% fix scale
%
% Revision 1.2  2007/04/05 21:17:05  arno
% fixed log freuqncies
%
% Revision 1.1  2007/01/26 18:10:15  arno
% Initial revision
%
% Revision 1.24  2006/11/23 00:26:22  arno
% cosmetic change
%
% Revision 1.23  2006/11/22 20:22:38  arno
% filter if necessary
%
% Revision 1.22  2006/11/22 19:34:22  arno
% empty plot
%
% Revision 1.21  2006/11/16 00:07:31  arno
% fix title for topoplot
%
% Revision 1.20  2006/11/15 21:58:49  arno
% plotting titles
%
% Revision 1.19  2006/11/09 23:55:18  arno
% fix last change
%
% Revision 1.18  2006/11/09 23:27:11  arno
% figure titles
%
% Revision 1.17  2006/11/09 23:21:21  arno
% add component name
%
% Revision 1.16  2006/11/09 22:04:39  arno
% ERSP plotting
%
% Revision 1.15  2006/11/03 03:01:17  arno
% allow plotting specific time-freq point
%
% Revision 1.14  2006/11/02 22:14:15  arno
% g. -> opt.
%
% Revision 1.13  2006/11/02 22:13:04  arno
% Limit to ERSP
%
% Revision 1.12  2006/11/02 21:59:18  arno
% input option
%
% Revision 1.11  2006/11/02 21:54:51  arno
% same
%
% Revision 1.10  2006/11/02 21:53:06  arno
% condstats -> condstat
%
% Revision 1.9  2006/10/10 23:50:54  scott
% replaced ?? defaults with defaults from finputcheck()
%
% Revision 1.8  2006/10/09 23:51:45  scott
% some more help edits
%
% Revision 1.7  2006/10/04 23:55:34  toby
% Bug fix courtesy Bas de Kruif
%
% Revision 1.6  2006/10/03 21:46:11  scott
% edit help msg -- ?? remain... -sm
%
% Revision 1.5  2006/10/03 16:29:27  scott
% edit
%
% Revision 1.4  2006/10/03 16:24:12  scott
% help message eidts.  ARNO - SEE ??s   -sm
%
% Revision 1.3  2006/10/02 11:41:13  arno
% wrote documentation
%

% Copyright (C) 2006 Arnaud Delorme
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

function [pgroup, pcond, pinter] = std_plottf(timevals, freqs, data, varargin)

pgroup = [];
pcond  = [];
pinter = [];
if nargin < 2
    help std_plottf;
    return;
end;

opt = finputcheck( varargin, { 'channels'       'cell'   []              {};
                               'titles'         'cell'   []              {};
                               'caxis'          'real'   []              [];
                               'ersplim'        'real'   []              []; % same as above
                               'itclim'         'real'   []              []; % same as above
                               'ylim'           'real'   []              [];
                               'condnames'      'cell'   []              {}; % just for titles
                               'groupnames'     'cell'   []              {}; % just for titles
                               'compinds'       'cell'   []              {};
                               'tftopoopt'      'cell'   []              {};
                               'threshold'      'real'   []              NaN;
                               'unitx'          'string' []              'ms'; % just for titles
                               'subject'        'string' []              '';   % just for titles
                               'chanlocs'       'struct' []              struct('labels', {});
                               'freqscale'      'string' { 'log' 'linear' 'auto' }  'auto';
                               'groupstats'     'cell'   []              {};
                               'condstats'      'cell'   []              {};
                               'interstats'     'cell'   []              {};                               
                               'maskdata'       'string' { 'on' 'off' }   'off';
                               'datatype'       'string' { 'ersp' 'itc' }    'ersp';
                               'plotmode'       'string' { 'normal' 'condensed' }  'normal' }, 'std_plottf');
if isstr(opt), error(opt); end;
if all(all(cellfun('size', data, 3)==1))               opt.singlesubject = 'on'; end;
if ~isempty(opt.compinds), if length(opt.compinds{1}) > 1, opt.compinds = {}; end; end;
if ~isempty(opt.groupstats) & ~isempty(opt.condstats) & strcmpi(opt.maskdata, 'on')
    disp('Cannot use ''maskdata'' option with both condition stat. and group stat. on');
    disp('Disabling statistics');
    opt.groupstats = {}; opt.condstats = {}; opt.maskdata = 'off'; 
end;
if ~isempty(opt.ersplim), opt.caxis = opt.ersplim; end;
if ~isempty(opt.itclim), opt.caxis = opt.itclim; end;
onecol  = { 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' 'b' };
manycol = { 'b' 'r' 'g' 'k' 'c' 'y' };

nc = size(data,1);
ng = size(data,2);
if nc >= ng, opt.transpose = 'on';
else         opt.transpose = 'off';
end;
if isempty(opt.condnames)
    for c=1:nc, opt.condnames{c} = sprintf('Cond. %d', c); end;
    if nc == 1, opt.condnames = { '' }; end;
end;
if isempty(opt.groupnames)
    for g=1:ng, opt.groupnames{g} = sprintf('Group. %d', g); end;
    if ng == 1, opt.groupnames = { '' }; end;
end;

% test log frequencies
% --------------------
if length(freqs) > 2 & strcmpi(opt.freqscale, 'auto')
    midfreq = (freqs(3)+freqs(1))/2;
    if midfreq*.9999 < freqs(2) & midfreq*1.0001 > freqs(2), opt.freqscale = 'linear';
    else                                                     opt.freqscale = 'log';
    end;
end;

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
            end;
        end;
    end;
    meanplot = meanplot/count;
    options = { 'chanlocs', opt.chanlocs, 'electrodes', 'off', 'cbar', 'on', ...
            'cmode', 'separate', opt.tftopoopt{:} };       
        
    if strcmpi(opt.freqscale, 'log'), options = { options{:} 'logfreq', 'native' }; end;
    tftopo( meanplot', timevals, freqs, 'title', [ 'Mean ' upper(opt.datatype) ' for all group/cond' ], options{:}); 
    currentHangle = gca;
    if ~isempty( opt.caxis )
        caxis( currentHangle, opt.caxis )
    end
    colorbarHandle = cbar;
    title(colorbarHandle,'dB');
    axes(currentHangle); 
    return; 
end;

% plotting paramters
% ------------------
if ng > 1 & ~isempty(opt.groupstats), addc = 1; else addc = 0; end;
if nc > 1 & ~isempty(opt.condstats  ), addr = 1; else addr = 0; end;

% compute significance mask
% --------------------------
if ~isempty(opt.interstats), pinter = opt.interstats{3}; end;

if ~isnan(opt.threshold) & ( ~isempty(opt.groupstats) | ~isempty(opt.condstats) )    
    pcondplot  = opt.condstats;
    pgroupplot = opt.groupstats;
    maxplot = 1;
else
    warning off;
    for ind = 1:length(opt.condstats),  pcondplot{ind}  = -log10(opt.condstats{ind}); end;
    for ind = 1:length(opt.groupstats), pgroupplot{ind} = -log10(opt.groupstats{ind}); end;
    if ~isempty(pinter), pinterplot = -log10(pinter); end;
    maxplot = 3;
    warning on;
end;

% -------------------------------
% masking for significance of not
% -------------------------------
statmask = 0;
if strcmpi(opt.maskdata, 'on') & ~isnan(opt.threshold) & ...
        (~isempty(opt.condstats) | ~isempty(opt.condstats))
    addc = 0; addr = 0; statmask = 1;
end;

% -------------------------
% plot time/frequency image
% -------------------------
options = { 'chanlocs', opt.chanlocs, 'electrodes', 'off', 'cbar', 'off', ...
            'cmode', 'separate', opt.tftopoopt{:} };
if strcmpi(opt.freqscale, 'log'), options = { options{:} 'logfreq', 'native' }; end;

figure('color', 'w');
tmpc = [inf -inf];
for c = 1:nc
    for g = 1:ng
        hdl(c,g) = mysubplot(nc+addr, ng+addc, g + (c-1)*(ng+addc), opt.transpose);
        if ~isempty(data{c,g})
            tmpplot = mean(data{c,g},3);
            if statmask, 
                if ~isempty(opt.condstats), tmpplot(find(pcondplot{g}(:) == 0)) = 0;
                else                        tmpplot(find(pgroupplot{c}(:) == 0)) = 0;
                end;
            end;

            tftopo( tmpplot', timevals, freqs, 'title', opt.titles{c,g}, options{:}); 
                
            if isempty(opt.caxis) & ~isempty(tmpc)
                warning off;
                tmpc = [ min(min(tmpplot(:)), tmpc(1)) max(max(tmpplot(:)), tmpc(2)) ];
                warning on;
            else 
                if ~isempty(opt.caxis)
                    caxis(opt.caxis);
                end;
            end;

            if c > 1
                ylabel(''); 
            end;
        end;
    
        % statistics accross groups
        % -------------------------
        if g == ng & ng > 1 & ~isempty(opt.groupstats) & ~statmask
            hdl(c,g+1) = mysubplot(nc+addr, ng+addc, g + 1 + (c-1)*(ng+addc), opt.transpose);
            tftopo( pgroupplot{c}', timevals, freqs, 'title', opt.titles{c,g+1}, options{:});
            caxis([-maxplot maxplot]);
        end;

    end;
end;


for g = 1:ng
    % statistics accross conditions
    % -----------------------------
    if ~isempty(opt.condstats) & ~statmask & nc > 1
        hdl(nc+1,g) = mysubplot(nc+addr, ng+addc, g + c*(ng+addc), opt.transpose);
        tftopo( pcondplot{g}', timevals, freqs, 'title', title(opt.titles{nc+1,g}), options{:});
        caxis([-maxplot maxplot]);
    end;
end;
            ylabel('');  % nima
% color scale
% -----------
if isempty(opt.caxis)
    tmpc = [-max(abs(tmpc)) max(abs(tmpc))];
    for c = 1:nc
        for g = 1:ng
            axes(hdl(c,g));
            if ~isempty(tmpc)
                caxis(tmpc);
            end;
        end;
    end;
end;

% statistics accross group and conditions
% ---------------------------------------
if ~isempty(opt.groupstats) & ~isempty(opt.condstats) & ng > 1 & nc > 1
    hdl(nc+1,ng+1) = mysubplot(nc+addr, ng+addc, g + 1 + c*(ng+addr), opt.transpose);
    tftopo( pinterplot',  timevals, freqs, 'title', opt.titles{nc+1,ng+1}, options{:});
    caxis([-maxplot maxplot]);
    ylabel('');
end;    

% color bars
% ----------
axes(hdl(nc,ng)); 
cbar_standard(opt.datatype, ng); 
if nc ~= size(hdl,1) | ng ~= size(hdl,2)
    axes(hdl(end,end));
    cbar_signif(ng, maxplot);
end;

% mysubplot (allow to transpose if necessary)
% -------------------------------------------
function hdl = mysubplot(nr,nc,ind,transp);

    r = ceil(ind/nc);
    c = ind -(r-1)*nc;
    if strcmpi(transp, 'on'), hdl = subplot(nc,nr,(c-1)*nr+r);
    else                      hdl = subplot(nr,nc,(r-1)*nc+c);
    end;

% colorbar for ERSP and scalp plot
% --------------------------------
function cbar_standard(datatype, ng);
    pos = get(gca, 'position');
    tmpc = caxis;
    fact = fastif(ng == 1, 40, 20);
    tmp = axes('position', [ pos(1)+pos(3)+max(pos(3)/fact,0.006) pos(2) max(pos(3)/fact,0.01) pos(4) ]);  
    set(gca, 'unit', 'normalized');
    if strcmpi(datatype, 'itc')
         cbar(tmp, 0, tmpc, 10); ylim([0.5 1]);
         title('ITC');
    else cbar(tmp, 0, tmpc, 5);title('dB');
    end;
    

% colorbar for significance
% -------------------------
function cbar_signif(ng, maxplot);
    pos = get(gca, 'position');
    tmpc = caxis;
    fact = fastif(ng == 1, 40, 20);
    tmp = axes('position', [ pos(1)+pos(3)+max(pos(3)/fact,0.006) pos(2) max(pos(3)/fact,0.01) pos(4) ]);  
    map = colormap;
    n = size(map,1);
    cols = [ceil(n/2):n]';
    image([0 1],linspace(0,maxplot,length(cols)),[cols cols]);
    %cbar(tmp, 0, tmpc, 5);
    tick = linspace(0, maxplot, maxplot+1);
    set(gca, 'ytickmode', 'manual', 'YAxisLocation', 'right', 'xtick', [], ...
        'ytick', tick, 'yticklabel', round(10.^-tick*1000)/1000);
    xlabel('');

% rapid filtering for ERP
% -----------------------
function tmpdata2 = myfilt(tmpdata, lowpass, highpass, factor, filtertype)

    tmpdata2 = reshape(tmpdata, size(tmpdata,1), size(tmpdata,2)*size(tmpdata,3)*size(tmpdata,4));
    tmpdata2 = eegfiltfft(tmpdata2',lowpass, highpass, factor, filtertype)';
    tmpdata2 = reshape(tmpdata2, size(tmpdata,1), size(tmpdata,2), size(tmpdata,3), size(tmpdata,4));
