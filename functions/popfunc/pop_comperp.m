% pop_comperp() - Compute the grand average ERP waveforms of multiple datasets,
%                 with optional ERP subtraction.
% Usage:
%       >> pop_comperp( ALLEEG, flag );      % pop-up window mode
%       >> [erp1 erp2 erpsub time sig] = pop_comperp( ALLEEG, flag, ...
%                                   datadd, datsub, 'key', 'val', ...);
% Inputs:
%   ALLEEG  - Array of EEG datasets
%   flag    - [0|1] (1) Use raw data or (0) ICA components. {default: 1}
%   datadd  - [integer array] List of datasets to sum to make an ERP 
%             grand average.
%   datsub  - [integer array] List of datasets to subtract to make an ERP
%             grand average. This option is to be used to compare
%             sub-conditions; each array must have the same number of 
%             elements as the respective 'datadd' datasets. 
%             e.g., The first 'datsub' dataset is subtracted from first 
%             'datadd' dataset. The two datasets may normally contain data 
%             from the same subject under different experiemtal conditions.
%
% Optional inputs:
%   'chans'    - [integer array] Vector of channel or component indices. 
%                {default: all}.
%   'title'    - [string] Plot title. {default: none}
%   'alpha'    - [float] Apply t-test for p=alpha (0<alpha<1). Use paired if
%                t-test datasub is not empty (two-tailed). If data is empty,
%                use t-test against a 0 mean dataset with same variance (two-
%                tailed). Significan time regions are highlighted in the
%                data plots.
%   'geom'     - ['scalp'|'array'] Plot erps in a scalp array (plottopo())
%                or as a rectangular array (plotdata()). Note: Components
%                cannot be plotted in a 'scalp' array.
%   'mode'     - ['ave'|'rms'] Plot grand average or RMS (root mean square)
%   'addavg'   - ['on'|'off'] Show average/RMS for datadd. Default: 'on' if
%                datsub empty, otherwise 'off'.
%   'addstd'   - ['on'|'off'] Show standard deviation for datadd. Default: 'on'
%                if datsub empty, otherwise 'off'.
%   'addall'   - ['on'|'off'] Show all erps for dataadd. Default: 'off'.
%   'subavg'   - ['on'|'off'] Show average/RMS for datsub. Default: 'off'.
%   'substd'   - ['on'|'off'] Show standard deviation for datsub. Default: 'off'.
%   'suball'   - ['on'|'off'] Show all erps for datasub. Default: 'off'.
%   'diffadd'  - ['on'|'off'] Show average/RMS for difference. Default: 'on'.
%   'diffstd'  - ['on'|'off'] Show standard deviation for difference.
%                Default is 'on'.
%   'diffall'  - ['on'|'off'] Show all erps for difference. Default: 'off'.
%   'tplotopt' - [cell array] 'key', val' plotting options for plottopo
%
% Output:
%   erp1   - Grand average (or rms) of 'datadd' datasets
%   erp2   - Grand average (or rms) of 'datsub' datasets
%   erpsub - Grand average (or rms) 'datadd' minus 'datsub' difference
%   times  - Array of epoch time indices
%   sig    - P significance values (chans,times). 
%
% Author: Arnaud Delorme, CNL / Salk Institute, 15 March 2003
%
% Note: t-test functions were adapted for matric preocessing from 
%       functions of C. Goutte. See description inside the code of this
%       function for more info.
%
% See also: eeglab(), plottopo()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 15 March 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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
% Revision 1.24  2003/12/04 17:49:06  arno
% same
%
% Revision 1.23  2003/12/04 17:47:26  arno
% eegfiltfft
%
% Revision 1.22  2003/12/04 15:34:27  arno
% fixing erp add
%
% Revision 1.21  2003/11/18 16:35:40  scott
% Enter ... -> ...
%
% Revision 1.20  2003/11/04 15:29:53  arno
% returning empty pvalues if alpha not set
%
% Revision 1.19  2003/07/29 19:10:31  arno
% same thing
%
% Revision 1.18  2003/07/29 19:01:50  arno
% handle empty ICA matrices
%
% Revision 1.17  2003/07/16 18:54:47  arno
% changing default
%
% Revision 1.16  2003/07/15 22:00:07  arno
% fixing backward compatibility problem
%
% Revision 1.15  2003/07/15 21:14:40  arno
% changed name of dataset
%
% Revision 1.14  2003/07/15 21:12:33  arno
% added backward compatibility
%
% Revision 1.13  2003/07/15 18:49:43  arno
% debug remaning problems
%
% Revision 1.12  2003/07/15 18:40:07  arno
% debugging difference ...
%
% Revision 1.11  2003/07/15 17:19:04  arno
% reprogrammed the function completelly to allow more options
%
% Revision 1.10  2003/05/10 17:59:21  arno
% debug command output
%
% Revision 1.9  2003/04/15 23:38:35  arno
% debuging last
%
% Revision 1.8  2003/04/15 17:10:14  arno
% adding low pass filtering capabilities
%
% Revision 1.7  2003/04/15 00:26:15  arno
% debug alpha
%
% Revision 1.6  2003/04/15 00:15:32  arno
% axis off
%
% Revision 1.5  2003/04/15 00:12:30  arno
% debuging for 1 dataset only
%
% Revision 1.4  2003/03/18 00:27:53  arno
% adding help button
%
% Revision 1.3  2003/03/17 23:42:31  arno
% adding lots of options
%
% Revision 1.2  2003/03/16 07:47:26  scott
% header msg
%
% Revision 1.1  2003/03/16 03:00:52  arno
% Initial revision
%

function [erp1, erp2, erpsub, time, pvalues] = pop_comperp( ALLEEG, flag, datadd, datsub, varargin);

erp1 = '';
if nargin < 1
   help pop_comperp;
   return;
end;
if isempty(ALLEEG)
    error('pop_comperp: cannot process empty sets of data');
end;
if nargin < 2
   flag = 1;
end;
allcolors = { 'b' 'r' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' ...
              'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm' 'r' 'b' 'g' 'c' 'm'};
erp1 = '';
if nargin < 3
    gtmp = [1.1 0.8 .18 .18 .18 0.1]; gtmp2 = [1.48 1.07 0.93];
    uigeom = { [2.6 0.95] gtmp gtmp gtmp [1] gtmp2 gtmp2 [1.48 0.21 1.79] gtmp2 [0.6 0.5 0.9] };
    commulcomp= ['if get(gcbo, ''value''),' ...
                 '    set(findobj(gcbf, ''tag'', ''multcomp''), ''enable'', ''on'');' ...
                 'else,' ...
                 '    set(findobj(gcbf, ''tag'', ''multcomp''), ''enable'', ''off'');' ...
                 'end;'];
	uilist = { { } ...
               { 'style' 'text' 'string' 'avg.        std.      all ERPs' } ...
               { 'style' 'text' 'string' 'Datasets to add (ex: 1 3 4):' } ...
               { 'style' 'edit' 'string' '' } ...
               { 'style' 'checkbox' 'string' '' 'value' 1 } ...
               { 'style' 'checkbox' 'string' '' } ...
               { 'style' 'checkbox' 'string' '' } { } ...
	           { 'style' 'text' 'string' 'Datasets to subtract (ex: 5 6 7):' } ...
               { 'style' 'edit' 'string' '' } ...
               { 'style' 'checkbox' 'string' '' 'value' 1 } ...
               { 'style' 'checkbox' 'string' '' } ...
               { 'style' 'checkbox' 'string' '' } { } ...
	           { 'style' 'text' 'string' 'Plot difference' } { } ...
               { 'style' 'checkbox' 'string' '' 'value' 1 } ...
               { 'style' 'checkbox' 'string' '' } ...
               { 'style' 'checkbox' 'string' '' } { } ...
               { } ...
	           { 'style' 'text' 'string' fastif(flag, 'Channels subset ([]=all):', ...
                                                  'Components subset ([]=all):') } ...
               { 'style' 'edit' 'string' '' } { } ...
	           { 'style' 'text' 'string' 'Highlight significant regions (.01 -> p=.01)' } ...
               { 'style' 'edit' 'string' '' } { } ...
	           { 'style' 'text' 'string' 'Use RMS instead of average (check):' } { 'style' 'checkbox' 'string' '' } { } ...
	           { 'style' 'text' 'string' 'Low pass (Hz) (for display only)' } ...
               { 'style' 'edit' 'string' '' } { } ...
               { 'style' 'text' 'string' 'Plottopo options (''key'', ''val''):' } ...
               { 'style' 'pushbutton' 'string' 'Help' 'callback', 'pophelp(''plottopo'')' } ...
               { 'style' 'edit' 'string' '''ydir'', -1' }};
    
    % remove geometry textbox for ICA components
    result = inputgui( uigeom, uilist, 'pophelp(''pop_comperp'')', 'ERP grand average/RMS - pop_comperp()');
    if length(result) == 0, return; end;

    %decode parameters list
    options = {};
    datadd = eval( [ '[' result{1} ']' ]);
    if result{2},  options = { options{:} 'addavg'  'on' }; else, options = { options{:} 'addavg'  'off' }; end;
    if result{3},  options = { options{:} 'addstd'  'on' }; else, options = { options{:} 'addstd'  'off' }; end;
    if result{4},  options = { options{:} 'addall'  'on' }; end;
    datsub = eval( [ '[' result{5} ']' ]);
    if result{6},  options = { options{:} 'subavg'  'on' }; end;
    if result{7},  options = { options{:} 'substd'  'on' }; end;
    if result{8},  options = { options{:} 'suball'  'on' }; end;
    if result{9},  options = { options{:} 'diffavg' 'on' }; else, options = { options{:} 'diffavg' 'off' }; end;
    if result{10}, options = { options{:} 'diffstd' 'on' }; else, options = { options{:} 'diffstd' 'off' }; end;
    if result{11}, options = { options{:} 'diffall' 'on' }; end;
    
    if result{12},           options = { options{:} 'chans' eval( [ '[' result{12} ']' ]) }; end;
    if ~isempty(result{13}), options = { options{:} 'alpha' str2num(result{13}) }; end;
    if result{14},           options = { options{:} 'mode' 'rms' }; end;
    if ~isempty(result{15}), options = { options{:} 'lowpass' str2num(result{15}) }; end;
    if ~isempty(result{16}), options = { options{:} 'tplotopt' eval([ '{ ' result{16} ' }' ]) }; end; 
else 
    options = varargin;
end;

% decode inputs
% -------------
if isempty(datadd), error('First edit box (datasets to add) can not be empty'); end;
g = finputcheck( options, ... 
                 { 'chans'    'integer'  [0:ALLEEG(datadd(1)).nbchan] 0;
                   'title'    'string'   []               '';
                   'alpha'    'float'    []               [];
                   'geom'     'string'  {'scalp' 'array'} fastif(flag, 'scalp', 'array');
                   'addstd'   'string'  {'on' 'off'}      fastif(isempty(datsub), 'on', 'off');
                   'substd'   'string'  {'on' 'off'}     'off';
                   'diffstd'  'string'  {'on' 'off'}     'on';
                   'addavg'   'string'  {'on' 'off'}     fastif(isempty(datsub), 'on', 'off');
                   'subavg'   'string'  {'on' 'off'}     'off';
                   'diffavg'  'string'  {'on' 'off'}     'on';
                   'addall'   'string'  {'on' 'off'}     'off';
                   'suball'   'string'  {'on' 'off'}     'off';
                   'diffall'  'string'  {'on' 'off'}     'off';
                   'std'      'string'  {'on' 'off' 'none'}     'none';
                   'diffonly' 'string'  {'on' 'off' 'none'}     'none';
                   'allerps'  'string'  {'on' 'off' 'none'}     'none';
                   'lowpass'  'float'    [0 Inf]         [];
                   'tplotopt' 'cell'     []              {};
                   'mode'     'string'  {'ave' 'rms'}    'ave';
                   'multcmp'  'integer'  [0 Inf]         [] });
if isstr(g), error(g); end;

figure;
try, icadefs; set(gcf, 'color', BACKCOLOR); axis off; catch, end;

% backward compatibility of param
% -------------------------------
if ~strcmpi(g.diffonly, 'none')
    if strcmpi(g.diffonly, 'off'), g.addavg = 'on'; g.subavg = 'on'; end;
end;
if ~strcmpi(g.allerps, 'none')
    if isempty(datsub)
         g.addall  = g.allerps;
    else g.diffall = g.allerps;
    end;
end;
if ~strcmpi(g.std, 'none')
    if isempty(datsub)
         g.addstd  = g.std;
    else g.diffstd = g.std;
    end;
end;

% check consistency
% -----------------
if length(datsub) > 0 & length(datadd) ~= length(datsub)
    error('The number of component to subtract must be the same as the number of components to add');
end;
regions = {};
pnts   = ALLEEG(datadd(1)).pnts;
srate  = ALLEEG(datadd(1)).srate;
xmin   = ALLEEG(datadd(1)).xmin;
xmax   = ALLEEG(datadd(1)).xmax;
nbchan = ALLEEG(datadd(1)).nbchan;
chanlocs = ALLEEG(datadd(1)).chanlocs;
for index = union(datadd, datsub)
    if ALLEEG(index).pnts ~= pnts,     error(['Dataset '  int2str(index) ' does not have the same number of points as others']); end;
    if ALLEEG(index).xmin ~= xmin,     error(['Dataset '  int2str(index) ' does not have the same xmin as others']); end;
    if ALLEEG(index).xmax ~= xmax,     error(['Dataset '  int2str(index) ' does not have the same xmax as others']); end;
    if ALLEEG(index).nbchan ~= nbchan, error(['Dataset '  int2str(index) ' does not have the same number of channels as others']); end;
end;
    
% compute ERPs for add
% --------------------
for index = 1:length(datadd)
    TMPEEG = eeg_checkset(ALLEEG(datadd(index)));
    if flag == 1, erp1ind(:,:,index)  = mean(TMPEEG.data,3);
    elseif isempty(TMPEEG.icaact)
        tmpica              =  reshape((TMPEEG.icaweights*TMPEEG.icasphere)*TMPEEG.data(:,:), ...
                                       size(TMPEEG.icaweights,1), size(TMPEEG.data,2), size(TMPEEG.data,3));
        erp1ind(:,:,index)  = mean(tmpica,3);
    else          
        tmpica              =  reshape((TMPEEG.icaweights*TMPEEG.icasphere)*TMPEEG.data(:,:), ...
                                       size(TMPEEG.icaweights,1), size(TMPEEG.data,2), size(TMPEEG.data,3));
        erp1ind(:,:,index)  = mean(tmpica,3);
        erp1ind(:,:,index)  = mean(TMPEEG.icaact,3);
    end;
    addnames{index} = [ '#' int2str(datadd(index)) ' ' TMPEEG.setname ' (n=' int2str(TMPEEG.trials) ')' ];
    clear TMPEEG;
end;

% optional: subtract
% ------------------
colors = {}; % color aspect for curves
allcolors = { 'b' 'r' 'g' 'c' 'm' 'y' [0 0.5 0] [0.5 0 0] [0 0 0.5] [0.5 0.5 0] [0 0.5 0.5] [0.5 0 0.5] [0.5 0.5 0.5] };
allcolors = { allcolors{:} allcolors{:} allcolors{:} allcolors{:} allcolors{:} allcolors{:} };
allcolors = { allcolors{:} allcolors{:} allcolors{:} allcolors{:} allcolors{:} allcolors{:} };
if length(datsub) > 0 % dataset to subtract

    % compute ERPs for sub
    % --------------------
    for index = 1:length(datsub)
        TMPEEG = eeg_checkset(ALLEEG(datsub(index)));
        if flag == 1, erp2ind(:,:,index)  = mean(TMPEEG.data,3);
        elseif isempty(TMPEEG.icaact)
        else          tmpica              =  reshape((TMPEEG.icaweights*TMPEEG.icasphere)*TMPEEG.data(:,:), ...
                                                     size(TMPEEG.icaweights,1), size(TMPEEG.data,2), size(TMPEEG.data,3));
                      erp2ind(:,:,index)  = mean(tmpica,3);
                      erp2ind(:,:,index)  = mean(TMPEEG.icaact,3);
        end;
        subnames{index} = [ '#' int2str(datsub(index)) ' ' TMPEEG.setname '(n=' int2str(TMPEEG.trials) ')' ];
        clear TMPEEG
    end;
    
    l1 = size(erp1ind,3);
    l2 = size(erp2ind,3);
    allcolors1 = allcolors(3:l1+2);
    allcolors2 = allcolors(l1+3:l1+l2+3);
    allcolors3 = allcolors(l1+l2+3:end);
    [erps1, erpstd1, colors1, colstd1, legend1] = preparedata( erp1ind        , g.addavg , g.addstd , g.addall , g.mode, 'Add ' , addnames, 'b', allcolors1 );
    [erps2, erpstd2, colors2, colstd2, legend2] = preparedata( erp2ind        , g.subavg , g.substd , g.suball , g.mode, 'Sub ' , subnames, 'r', allcolors2 );
    [erps3, erpstd3, colors3, colstd3, legend3] = preparedata( erp1ind-erp2ind, g.diffavg, g.diffstd, g.diffall, g.mode, 'Diff ', ...
                                                      { addnames subnames }, 'k', allcolors3 );
    
    % handle special case of std
    % --------------------------
    erptoplot  = [ erps1 erps2 erps3 erpstd1 erpstd2 erpstd3 ];
    colors     = { colors1{:} colors2{:} colors3{:} colstd1{:} colstd2{:} colstd3{:}};
    legend     = { legend1{:} legend2{:} legend3{:} };
    
    % highlight significant regions
    % -----------------------------
    if ~isempty(g.alpha)
        pvalues = pttest(erp1ind, erp2ind, 3);
        regions = p2regions(pvalues, g.alpha, [xmin xmax]*1000);
    else 
        pvalues= [];
    end;
    
else
    [erptoplot, erpstd, colors, colstd, legend] = preparedata( erp1ind, g.addavg, g.addstd, g.addall, g.mode, '', addnames, 'k', allcolors);
    erptoplot = [ erptoplot erpstd ];
    colors    = { colors{:} colstd{:} };
    
    % highlight significant regions
    % -----------------------------
    if ~isempty(g.alpha)
        pvalues = ttest(erp1ind, 0, 3);
        regions = p2regions(pvalues, g.alpha, [xmin xmax]*1000);
    else 
        pvalues= [];
    end;
end;
    
% lowpass data
% ------------
if ~isempty(g.lowpass)
    if exist('filtfilt') == 2
        erptoplot = eegfilt(erptoplot, srate, 0, g.lowpass);
    else
        erptoplot = eegfiltfft(erptoplot, srate, 0, g.lowpass);
    end;
end;
if strcmpi(g.geom, 'array') | flag == 0, chanlocs = []; end;

% plot data
% ---------
plottopo( erptoplot, 'chanlocs', chanlocs, 'frames', pnts, ...
          'limits', [xmin xmax 0 0]*1000, 'title', g.title, 'colors', colors, ...
          'chans', g.chans, 'legend', legend, 'regions', regions, g.tplotopt{:});
times = linspace(xmin, xmax, pnts);

if nargin < 3 & nargout == 1
    erp1 = sprintf('pop_comperp( %s, %d, %s);', inputname(1), ...
                  flag, vararg2str({ datadd datsub options{:} }) );
end;
return;

% convert significance values to alpha
% ------------------------------------
function regions = p2regions( pvalues, alpha, limits);
    
    for index = 1:size(pvalues,1)
        signif = diff([1 pvalues(index,:) 1] < alpha);
        pos    = find([signif] > 0);
        pos    = pos/length(pvalues)*(limits(2) - limits(1))+limits(1);
        neg    = find([signif(2:end)] < 0);
        neg    = neg/length(pvalues)*(limits(2) - limits(1))+limits(1);
        if length(pos) ~= length(neg), signif, pos, neg, error('Region error'); end;
        regions{index} = [neg;pos];
    end;
    
% process data
% ------------
function [erptoplot, erpstd, colors, colstd, legend] = preparedata( erpind, plotavg, plotstd, plotall, mode, tag, dataset, coloravg, allcolors);

    colors    = {};
    legend    = {};
    erptoplot = [];
    erpstd    = [];
    colstd    = {};

    % plot individual differences
    % ---------------------------
    if strcmpi(plotall, 'on')
        erptoplot = [ erptoplot erpind(:,:) ];
        for index=1:size(erpind,3)
            if iscell(dataset)
                if strcmpi(tag, 'Diff ')
                    legend = { legend{:} [ dataset{1}{index} ' - ' dataset{2}{index} ] };
                else
                    legend = { legend{:} dataset{index} };
                end;
            else
                legend = { legend{:} [ 'Dataset ' int2str(dataset(index)) ] };
            end;
            colors = { colors{:}  allcolors{index} };
        end;
    end;
    
    % plot average
    % ------------
    if strcmpi( plotavg, 'on')
        if strcmpi(mode, 'ave')
             granderp    = mean(erpind,3);
             legend      = { legend{:} [ tag 'Average' ] };
        else granderp    = sqrt(mean(erpind.^2,3));
             legend      = { legend{:} [ tag 'RMS' ] };
        end;
        colors    = { colors{:}  {coloravg 'linewidth' 2 }};
        erptoplot = [ erptoplot granderp];
    end;

    % plot standard deviation
    % -----------------------
    if strcmpi(plotstd, 'on')
        if strcmpi(plotavg, 'on')
            std1      = std(erpind, [], 3);
            erptoplot = [ erptoplot granderp+std1 ];
            erpstd    = granderp-std1;
            legend    = { legend{:} [ tag 'Standard dev.' ] };
            colors    = { colors{:} { coloravg 'linestyle' ':' } };
            colstd    = { { coloravg 'linestyle' ':' } };
        else 
            disp('Warning: cannot show standard deviation without showing average');
        end;
    end;

% ------------------------------------------------------------------
    
function [p, t, df] = pttest(d1, d2, dim)
%PTTEST Student's paired t-test.
%       PTTEST(X1, X2) gives the probability that Student's t
%       calculated on paired data X1 and X2 is higher than
%       observed, i.e. the "significance" level. This is used
%       to test whether two paired samples have significantly
%       different means.
%       [P, T] = PTTEST(X1, X2) gives this probability P and the
%       value of Student's t in T. The smaller P is, the more
%       significant the difference between the means.
%       E.g. if P = 0.05 or 0.01, it is very likely that the
%       two sets are sampled from distributions with different
%       means.
%
%       This works for PAIRED SAMPLES, i.e. when elements of X1
%       and X2 correspond one-on-one somehow.
%       E.g. residuals of two models on the same data.
%       Ref: Press et al. 1992. Numerical recipes in C. 14.2, Cambridge.

if size(d1,dim) ~= size(d2, dim)
   error('PTTEST: paired samples must have the same number of elements !')
end
if size(d1,dim) == 1
    close; error('Cannot compute paired t-test for a single ERP difference')
end; 
a1 = mean(d1, dim);
a2 = mean(d2, dim);
v1 = std(d1, [], dim).^2;
v2 = std(d2, [], dim).^2;
n1 = size(d1,dim);
df = n1 - 1;
disp(['Computing t-values, df:' int2str(df) ]);

d1 = d1-repmat(a1, [ones(1,dim-1) size(d1,3)]);
d2 = d2-repmat(a2, [ones(1,dim-1) size(d2,3)]);
%cab = (x1 - a1)' * (x2 - a2) / (n1 - 1);
cab = mean(d1.*d2,3)/(n1-1);
% use abs to avoid numerical errors for very similar data
% for which v1+v2-2cab may be close to 0.
t = (a1 - a2) ./ sqrt(abs(v1 + v2 - 2 * cab) / n1) ;
p = betainc( df ./ (df + t.*t), df/2, 0.5) ;

% ------------------------------------------------------------------

function [p, t] = ttest(d1, d2, dim)
%TTEST Student's t-test for equal variances.
%       TTEST(X1, X2) gives the probability that Student's t
%       calculated on data X1 and X2, sampled from distributions
%       with the same variance, is higher than observed, i.e.
%       the "significance" level. This is used to test whether
%       two sample have significantly different means.
%       [P, T] = TTEST(X1, X2) gives this probability P and the
%       value of Student's t in T. The smaller P is, the more
%       significant the difference between the means.
%       E.g. if P = 0.05 or 0.01, it is very likely that the
%       two sets are sampled from distributions with different
%       means.
%
%       This works if the samples are drawn from distributions with
%       the SAME VARIANCE. Otherwise, use UTTEST.
%
%See also: UTTEST, PTTEST.
if size(d1,dim) == 1
    close; error('Cannot compute t-test for a single ERP')
end; 
a1 = mean(d1, dim);
v1 = std(d1, [], dim).^2;
n1 = size(d1,dim);
if length(d2) == 1 & d2 == 0
    a2 = 0;
    n2 = n1;
    df = n1 + n2 - 2;
    pvar = (2*(n1 - 1) * v1) / df ;
else 
    a2 = mean(d2, dim);
    v2 = std(d2, [], dim).^2;
    n2 = size(d2,dim);
    df = n1 + n2 - 2;
    pvar = ((n1 - 1) * v1 + (n2 - 1) * v2) / df ;     
end;
disp(['Computing t-values, df:' int2str(df) ]);

t = (a1 - a2) ./ sqrt( pvar * (1/n1 + 1/n2)) ;
p = betainc( df ./ (df + t.*t), df/2, 0.5) ;


