% pop_corrmap - calls corrmap() function which compute correlations between IC template from a single
%               dataset and ICs from different datasets, clusters ICs above correlation threshold and 
%               displays summary plot containing topographical maps for
%               clustered ICs, average map, correlation distribution and 
%
% Usage:
%          >> pop_corrmap(STUDY, ALLEEG, n_tmp, index, 'key', 'val')
%
% Inputs:
%  STUDY - input STUDY structure
%  ALLEEG - input ALLEEG structure
%  n_tmp - number of the template dataset
%  index - index for component of template dataset that is going to be
%          correlated with all ICs from all datasets stored in ALLEEG
%
% Optional inputs:
%   'key','val' - optional corrmap() arguments (see >> help corrmap)
%
% See also:  corrmap(), corrmap_plot_v1(), corrmap_plot_v2()
%
% Author: F. Campos Viola, MRC-IHR, 20/07/2007

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) F. Campos Viola, MRC-IHR
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

function [CORRMAP,STUDY,ALLEEG,com] = pop_corrmap(STUDY,ALLEEG,n_tmp,index,varargin);

CORRMAP=[];
com = 'babbb'; % this initialization ensure that the function will return something
% if the user press the cancel button

if nargin < 2
    disp('No input parameters')
    disp('see >> help pop_corrmap.');
    return;
end;

if nargin < 3
    promptstr    = { 'Cluster name (default: not stored):' ...
        'Template Dataset (Select one):', ...
        'Template IC (Select one):', ...
        'Correlation threshold (default: auto):',...
        'Max number of ICs per dataset (default: 2):',...
        'Create EEG.badcomps (default: not stored):',...
        };

    inistr       = { '', 1, ...
        1,...
        'auto',...
        2,...
        'no'};
    
    result = inputdlg2( promptstr, 'Correlation between IC maps -- pop_corrmap()', 1, inistr, 'pop_corrmap');

    if size(result,1) == 0
      
    else
        [ chanlocs warn ] = eeg_mergelocs(ALLEEG.chanlocs);
        if warn
            error( [ 'Different channel montage order in some datasets.' 10 ...
                     'CORRMAP will only work with identical channel' 10 ...
                     'montage order.' ]);
        end;
        
        clname=result{1};
        n_tmp   = eval(  result{2} );
        index    = eval( result{3} );
        th=result{4};
        ics=eval( result{5} );
        badcomps=result{6};
        title=result{1};
        if ~isempty(title), title = [ 'Cluster ' title ]; end;
        resetclusters = 'off';
        
        sameicas  = std_findsameica(ALLEEG);
        datinds   = cellfun(@(x)(x(1)), sameicas);
        totalicas = 0;
        for datind = datinds, totalicas = totalicas + size(ALLEEG(datind).icaweights,1); end;
        
        if totalicas ~= length(STUDY.cluster(1).comps)
            ButtonName = questdlg2([ 'The CORRMAP plugin use by default all ICA' 10 ...
                'components (even if you have selected' 10 ...
                'them by residual variance). This requires' 10 ...
                'reinitialization of the component selection.' 10 ...
                'Alternatively, this function may return only' 10 ...
                'the components you selected by residual variance,' 10 ...
                'by pruning the results returned by CORRMAP.'], '', 'Cancel', 'Prune','Reinit','Prune');
            if strcmpi(ButtonName, 'Cancel'), return; end;
            if strcmpi(ButtonName, 'Reinit'), resetclusters = 'on'; end;
        end;
       
        [CORRMAP,STUDY,ALLEEG] = corrmap(STUDY,ALLEEG,n_tmp,index,'chanlocs',chanlocs,'th',th,'ics',ics,'title',title,'clname',clname,'badcomps',badcomps, 'resetclusters', resetclusters);
        % com = sprintf('pop_corrmap(STUDY,ALLEEG,%g,%g,''chanlocs'',''%s'',''th'',''%s'',''ics'',%g,''title'',''%s'',''clname'',''%s'',''badcomps'',''%s'');', n_tmp,index,chanlocs,th,ics,title,clname,badcomps);
        com = sprintf('pop_corrmap(STUDY,ALLEEG,%g, %g,''chanlocs'','''',''th'',''%s'',''ics'',%g,''title'',''%s'',''clname'',''%s'',''badcomps'',''%s'', ''resetclusters'',''%s'');', n_tmp,index,th,ics,title,clname,badcomps,resetclusters);
    end

else
    [CORRMAP,STUDY,ALLEEG] = corrmap(STUDY,ALLEEG,n_tmp,index,'plot', 'off', varargin{:});
end

