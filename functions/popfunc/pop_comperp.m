% pop_comperp() - Compute the grand average ERP waveforms of multiple datasets,
%                 with optional ERP subtraction.
% Usage:
%       >> pop_comperp( ALLEEG, flag );      % pop-up window mode
%       >> [erp1 erp2 erpsub time sig] = pop_comperp( ALLEEG, flag, ...
%                                   datadd, datsub, chansubset, title);
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
%   'std'      - ['on'|'off'] Show standard deviation. Default: 'off'.
%   'diffonly' - ['on'|'off'] When subtracting datasets, ('on') do not plot 
%                or ('off') do plot the ERP grand averages.
%   'allerps'  - ['on'|'off'] Show all erps. {default: 'off'}
%   'mode'     - ['ave'|'rms'] Plot grand average or RMS (root mean square)
%   'tplotopt' - [cell array] 'key', val' plotting options for topoplot
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
    checkboxgeom = [1.6 0.15 0.5];
    uigeom = { [ 2 1] [2 1] [2 1] [2 1] checkboxgeom checkboxgeom checkboxgeom checkboxgeom [1.2 0.5 1] };
    commulcomp= ['if get(gcbo, ''value''),' ...
                 '    set(findobj(gcbf, ''tag'', ''multcomp''), ''enable'', ''on'');' ...
                 'else,' ...
                 '    set(findobj(gcbf, ''tag'', ''multcomp''), ''enable'', ''off'');' ...
                 'end;'];
	uilist = { { 'style' 'text' 'string' 'Enter a list of datasets to add (ex: 1 3 4):' } ...
               { 'style' 'edit' 'string' '' } ...
	           { 'style' 'text' 'string' 'Enter a list of datasets to subtract (ex: 5 6 7):' } ...
               { 'style' 'edit' 'string' '' } ...
	           { 'style' 'text' 'string' fastif(flag, 'Enter a subset of channels ([]=all):', ...
                                                  'Enter a subset of components ([]=all):') } ...
               { 'style' 'edit' 'string' '' } ...
	           { 'style' 'text' 'string' 'p value to highligh significant regions (i.e. 0.01)' } ...
               { 'style' 'edit' 'string' '' } ...
	           { 'style' 'text' 'string' 'Show standard deviations:' } ...
               { 'style' 'checkbox' 'string' '' } { } ...
	           { 'style' 'text' 'string' 'Show all ERPs:' } ...
               { 'style' 'checkbox' 'string' '' } { } ...
	           { 'style' 'text' 'string' 'Show difference plot(s) only:' } ...
               { 'style' 'checkbox' 'string' '' } { } ...
	           { 'style' 'text' 'string' 'Check checkbox to use RMS instead of average:' } ...
               { 'style' 'checkbox' 'string' '' } { } ...
               { 'style' 'text' 'string' 'Topoplot options (''key'', ''val''):' } ...
               { 'style' 'pushbutton' 'string' 'Help' 'callback', 'pophelp(''plottopo'')' } ...
               { 'style' 'edit' 'string' '' }};
    
    % remove geometry textbox for ICA components
    result = inputgui( uigeom, uilist, 'pophelp(''pop_comperp'')', 'ERP grand average/RMS - pop_comperp()');
    if length(result) == 0, return; end;

    %decode parameters
    datadd = eval( [ '[' result{1} ']' ]);
    datsub = eval( [ '[' result{2} ']' ]);
    options = {};
    if result{3}, options = { options{:} 'chans' eval( [ '[' result{3} ']' ]) }; end;
    if ~isempty(result{4}), options = { options{:} 'alpha' result{4} }; end;
    if result{5}, options = { options{:} 'std' 'on' }; end;
    if result{6}, options = { options{:} 'allerps' 'on' }; end;
    if result{7}, options = { options{:} 'diffonly' 'on' }; end;
    if result{8}, options = { options{:} 'mode' 'rms' }; end;
    if ~isempty(result{9}), options = { options{:} 'tplotopt' eval([ '{ ' result{9} ' }' ]) }; end; 
else 
    options = varargin;
end;

% decode inputs
% -------------
if isempty(datadd), error('First edit box (datasets to add) can not be empty'); end;
g = finputcheck( options, ... 
                 { 'chans'    'integer'  [1:ALLEEG(datadd(1)).nbchan] [];
                   'title'    'string'   []               '';
                   'alpha'    'float'    []               [];
                   'geom'     'string'  {'scalp' 'array'} fastif(flag, 'scalp', 'array');
                   'std'      'string'  {'on' 'off'}     'off';
                   'diffonly' 'string'  {'on' 'off'}     'off';
                   'allerps'  'string'  {'on' 'off'}     'off';
                   'tplotopt' 'cell'     []              {};
                   'mode'     'string'  {'ave' 'rms'}    'ave';
                   'multcmp'  'integer'  [0 Inf]         [] });
if isstr(g), error(g); end;

figure;
try, icadefs; set(gcf, 'color', BACKCOLOR); axis off; catch, end;

% check consistency
% -----------------
if length(datsub) > 0 & length(datadd) ~= length(datsub)
    error('The number of component to subtract must be the same as the number of components to add');
end;
regions = {};
pnts = ALLEEG(datadd(1)).pnts;
xmin = ALLEEG(datadd(1)).xmin;
xmax = ALLEEG(datadd(1)).xmax;
nbchan = ALLEEG(datadd(1)).nbchan;
chanlocs = ALLEEG(datadd(1)).chanlocs;
for index = union(datadd, datsub)
    if ALLEEG(index).pnts ~= pnts,     error(['Dataset '  int2str(index) ' has not the same number of point as others']); end;
    if ALLEEG(index).xmin ~= xmin,     error(['Dataset '  int2str(index) ' has not the same xmin as others']); end;
    if ALLEEG(index).xmax ~= xmax,     error(['Dataset '  int2str(index) ' has not the same xmax as others']); end;
    if ALLEEG(index).nbchan ~= nbchan, error(['Dataset '  int2str(index) ' has not the same number of channels as others']); end;
end;
    
% compute grand average
% ---------------------
for index = 1:length(datadd)
    TMPEEG = eeg_checkset(ALLEEG(datadd(index)));
    if flag == 1, erp1ind(:,:,index)  = mean(TMPEEG.data,3);
    else          erp1ind(:,:,index)  = mean(TMPEEG.icaact,3);
    end;
    clear TMPEEG;
end;

% optional: subtract
% ------------------
colors = {}; % color aspect for curves
if length(datsub) > 0
    for index = 1:length(datsub)
        TMPEEG = eeg_checkset(ALLEEG(datsub(index)));
        if flag == 1, erp2ind(:,:,index)  = mean(TMPEEG.data,3);
        else          erp2ind(:,:,index)  = mean(TMPEEG.icaact,3);
        end;
        clear TMPEEG
    end;
    colors    = {{'k' 'linewidth' 2 }};
    if strcmpi(g.mode, 'ave')
         erpsub = mean(erp1ind-erp2ind,3);
         erp1   = mean(erp1ind,3);
         erp2   = mean(erp2ind,3);
         legend = { 'Avg. difference' };
    else erpsub = sqrt(mean((erp1ind-erp2ind).^2,3));
         erp1   = sqrt(mean(erp1ind.^2,3));
         erp2   = sqrt(mean(erp2ind.^2,3));
         legend = { 'RMS difference' };
    end;
    if ~isempty(g.alpha)
        pvalues = pttest(erp1ind, erp2ind, 3);
        regions = p2regions(pvalues, g.alpha, [xmin xmax]*1000);
    end;
    if strcmpi(g.diffonly, 'off') & strcmpi(g.allerps, 'on')
        error('can not plot individual averages with grand erps, set ''diffonly'' to ''on''');
    end;
    
    % plot individual differences
    % ---------------------------
    if strcmpi(g.allerps, 'on')
        erptoplot = [ erpsub erp1ind(:,:)-erp2ind(:,:) ];
        for index=1:size(erp1ind,3)
            legend{index+1} = [ 'Datasets ' int2str(datadd(index)) '-' int2str(datsub(index)) ];
        end;
    else 
        erptoplot = erpsub;   
    end;
    
    % plot grand averages
    % -------------------
    if strcmpi(g.diffonly, 'off')
        erptoplot = [ erpsub erp1 erp2];
        legend    = { legend{:} 'Dataset(s) added' 'Dataset(s) subtracted' };
        colors    = { colors{:} 'r' 'b' };
        if strcmpi(g.std, 'on')
            legend    = { legend{:} 'Std. diff.' 'Std. added' 'Std. subtracted' };
            stdsub    = std(erp1ind-erp2ind, [], 3);
            std1      = std(erp1ind, [], 3);
            std2      = std(erp2ind, [], 3);
            erptoplot = [ erptoplot erpsub+stdsub erp1+std1 erp2+std2 erpsub-stdsub erp1-std1 erp2-std2 ];
            colors    = { colors{:} 'k:' 'r:' 'b:' 'k:' 'r:' 'b:' };
        end;
    else 
        if strcmpi(g.std, 'on')
            stdsub    = std(erp1ind-erp2ind, [], 3);
            erptoplot = [ erptoplot erpsub+stdsub erpsub-stdsub ];
            legend    = { legend{:} 'Std. diff.' };
            colors    = { colors{:} 'k:' 'k:' };
        end;
    end;
else
    erpsub = []; erp2 = [];
    std1 = std(erp1ind, [], 3);
    if strcmpi(g.mode, 'ave')
         erp1   = mean(erp1ind,3);
         legend = { 'Avg.' };
    else erp1 = sqrt(mean(erp1ind.^2,3));
         legend = { 'RMS' };
    end;
    erptoplot = erp1;
    colors    = {{'k' 'linewidth' 2 }};

    % highlight significant regions
    % -----------------------------
    if ~isempty(g.alpha)
        pvalues = ttest(erp1ind, 0, 3);
        regions = p2regions(pvalues, g.alpha, [xmin xmax]*1000);
    end;

    % plot individual differences
    % ---------------------------
    if strcmpi(g.allerps, 'on')
        erptoplot = [ erptoplot erp1ind(:,:) ];
        for index=1:size(erp1ind,3)
            legend{index+1} = [ 'Dataset ' int2str(datadd(index)) ];
            colors    = { colors{:} allcolors{index} };
        end;
    end;
    
    % plot standard deviation
    % -----------------------
    if strcmpi(g.std, 'on')
        erptoplot = [ erptoplot erp1+std1 erp1-std1 ];
        legend    = { legend{:} 'Std.' };
        colors    = { colors{:} 'k:' 'k:' };
    end;

end;
    
if strcmpi(g.geom, 'array') | flag == 0, chanlocs = []; end;

plottopo( erptoplot, 'chanlocs', chanlocs, 'frames', pnts, ...
          'limits', [xmin xmax 0 0]*1000, 'title', g.title, 'colors', ...
          colors, 'legend', legend, 'regions', regions, g.tplotopt{:});
times = linspace(xmin, xmax, pnts);

if nargin < 3 & nargout == 1
    erp1 = sprintf('pop_comperp( %s, %d, [%s], [%s] %s);', inputname(1), ...
                  flag, num2str(datadd), num2str(datsub), vararg2str(options) );
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
disp(['Computing t-values, df:' int2str(df) ]);
a1 = mean(d1, dim);
a2 = mean(d2, dim);
v1 = std(d1, [], dim).^2;
v2 = std(d2, [], dim).^2;
n1 = size(d1,dim);
df = n1 - 1;

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


