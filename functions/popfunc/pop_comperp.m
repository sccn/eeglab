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
%   'alpha'    - [float] Apply t-test for p=alpha. If [float N], only mark 
%                for significance when N consecutive points are significant.
%                {default: []}
%   'geom'     - ['scalp'|'array'] Plot erps in a scalp array (plottopo())
%                or as a rectangular array (plotdata()). Note: Components
%                cannot be plotted in a 'scalp' array.
%   'std'      - ['on'|'off'] Show standard deviation. Default: 'off'.
%   'diffonly' - ['on'|'off'] When subtracting datasets, ('on') do not plot 
%                or ('off') do plot the ERP grand averages.
%   'allerps'  - ['on'|'off'] Show all erps. {default: 'off'}
%   'mode'     - ['ave'|'rms'] Plot grand average or RMS (root mean square)
%   'multcmp'  - [integer] Correct for multiple comparisons. Enter the number
%                of degrees of freedom (divides alpha by this number).
% Output:
%   erp1   - Grand average (or rms) of 'datadd' datasets
%   erp2   - Grand average (or rms) of 'datsub' datasets
%   erpsub - Grand average (or rms) 'datadd' minus 'datsub' difference
%   times  - Array of epoch time indices
%   sig    - Significance indicator (array of 0s and 1s of the same size as
%            erp1 or erpsub).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 15 March 2003
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
% Revision 1.1  2003/03/16 03:00:52  arno
% Initial revision
%

function [erp1, erp2, erpsub, time, sig] = pop_comperp( ALLEEG, flag, chanadd, chansub, varargin);

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

if nargin < 3
    checkboxgeom = [1.6 0.15 0.5];
    uigeom = { [ 2 1] [2 1] [2 1] [2 1] checkboxgeom checkboxgeom checkboxgeom checkboxgeom checkboxgeom checkboxgeom };
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
	           { 'style' 'text' 'string' 'Alpha for student two-tailed t-test ([]=none):' } ...
               { 'style' 'edit' 'string' '' } ...
	           { 'style' 'text' 'string' 'Correct for multiple comparisons:' } ...
               { 'style' 'checkbox' 'string' '' 'callback' commulcomp } ...
               { 'style' 'edit' 'string' num2str(ALLEEG(end).pnts) 'enable' 'off' 'tag' 'multcomp' } ...
	           { 'style' 'text' 'string' 'Show standard deviations:' } ...
               { 'style' 'checkbox' 'string' '' } { } ...
	           { 'style' 'text' 'string' 'Show all ERPs:' } ...
               { 'style' 'checkbox' 'string' '' } { } ...
	           { 'style' 'text' 'string' 'Show difference plot(s) only:' } ...
               { 'style' 'checkbox' 'string' '' } { } ...
	           { 'style' 'text' 'string' 'Check checkbox to use RMS instead of average:' } ...
               { 'style' 'checkbox' 'string' '' } { } ...
	           { 'style' 'text' 'string' 'Use rectangular geometry for plotting channels:' } ...
               { 'style' 'checkbox' 'string' '' } { } };
    
    % remove geometry textbox for ICA components
    if ~flag, uigeom(end) = []; uilist(end-2:end) = []; end;
    result = inputgui( uigeom, uilist, 'pophelp(''pop_comperp'')', 'ERP grand average/RMS - pop_comperp()');
    if length(result) == 0, return; end;

    %decode parameters
    datadd = eval( [ '[' result{1} ']' ]);
    datsub = eval( [ '[' result{2} ']' ]);
    options = {};
    if result{3}, options = { options{:} 'chans' eval( [ '[' result{3} ']' ]) }; end;
    if result{4}, options = { options{:} 'alpha' eval( [ '[' result{4} ']' ]) }; end;
    if result{5}, options = { options{:} 'multcmp' eval( result{6} ) }; end;
    if result{7}, options = { options{:} 'std' 'on' }; end;
    if result{8}, options = { options{:} 'allerps' 'on' }; end;
    if result{9}, options = { options{:} 'diffonly' 'on' }; end;
    if result{10}, options = { options{:} 'mode' 'rms' }; end;
    if length(result) == 11 & result{11}, options = { options{:} 'geom' 'array' }; end; 
else 
    options = varargin;
end;

% decode inputs
% -------------
if isempty(datadd), error('First edit box (datasets to add) can not be empty'); end;
g = finputcheck( options, ... 
                 { 'chans'    'integer'  [1:ALLEEG(datadd(1)).nbchan] [];
                   'title'    'string'   []               '';
                   'alpha'    'real'     [0 1]            [];
                   'geom'     'string'  {'scalp' 'array'} fastif(flag, 'scalp', 'array');
                   'std'      'string'  {'on' 'off'}     'off';
                   'diffonly' 'string'  {'on' 'off'}     'off';
                   'allerps'  'string'  {'on' 'off'}     'off';
                   'mode'     'string'  {'ave' 'rms'}    'ave';
                   'multcmp'  'integer'  [0 Inf]         [] });
if isstr(g), error(g); end;

figure;
try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end;

% check consistency
% -----------------
if length(datsub) > 0 & length(datadd) ~= length(datsub)
    error('The number of component to subtract must be the same as the number of components to add');
end;
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
if length(datsub) > 1
    for index = 1:length(datsub)
        TMPEEG = eeg_checkset(ALLEEG(datsub(index)));
        if flag == 1, erp2ind(:,:,index)  = mean(TMPEEG.data,3);
        else          erp2ind(:,:,index)  = mean(TMPEEG.icaact,3);
        end;
        clear TMPEEG
    end;
    std1 = std(erp1-erp2, [], 3);
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
    erptoplot = erpsub;
    if strcmpi(g.diffonly, 'off') & strcmpi(g.allerps, 'on')
        error('can not plot individual averages with grand erps, set ''diffonly'' to ''on''');
    end;
    if strcmpi(g.diffonly, 'off')
        erptoplot = [ erptoplot erp1 erp2];
        legend    = { legend{:} 'Dataset(s) added' 'Dataset(s) subtracted' };
    end;
    if strcmpi(g.allerps, 'on')
        erptoplot = [ erptoplot erp1ind(:,:)-erp2ind(:,:) ];
        for index=1:size(erp1ind,3)
            legend{index+1} = [ 'Datasets ' int2str(datadd(index)) '-' int2str(datsub(index)) ];
        end;
        colors    = {{'k' 'linewidth' 2 }};
    end;
else
    std1 = std(erp1ind, [], 3);
    if strcmpi(g.mode, 'ave')
         erp1   = mean(erp1ind,3);
         legend = { 'Avg.' };
    else erp1 = sqrt(mean(erp1ind.^2,3));
         legend = { 'RMS' };
    end;
    erptoplot = erp1;
    if strcmpi(g.allerps, 'on')
        erptoplot = [ erptoplot erp1ind(:,:) ];
        for index=1:size(erp1ind,3)
            legend{index+1} = [ 'Dataset ' int2str(datadd(index)) ];
        end;
        colors    = {{'k' 'linewidth' 2 }};
    end;
end;
    
if strcmpi(g.geom, 'array') | flag == 0, chanlocs = []; end;

plottopo( erptoplot, 'chanlocs', chanlocs, 'frames', pnts, ...
          'limits', [xmin xmax 0 0]*1000, 'title', g.title, 'colors', colors, 'legend', legend);

%plotdata( erptoplot, pnts, [xmin*1000 xmax*1000 0 0], g.title, g.chans );

%com = sprintf('figure; pop_comperp( %s, [%s], [%s], ''%s'');', inputname(1), num2str(setlist), num2str(chansubset), plottitle);
return;
