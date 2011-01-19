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
    promptstr    = { 'Template Dataset (Select one):', ...
        'Template IC (Select one):', ...
        'Correlation threshold (default: auto):',...
        'Number of ICs (default: 2):',...
        'Create EEG.badcomps (default: not stored):',...
        'Plot Title (default: no title):',...
        'Cluster name (default: not stored):'};

    inistr       = { 1, ...
        1,...
        'auto',...
        2,...
        'no',...
        '',...
        ''};
    
    result = inputdlg2( promptstr, 'Correlation between IC maps -- pop_corrmap()', 1, inistr, 'pop_corrmap');

    if size(result,1) == 0
      
    else
        chanlocs=eeg_mergelocs(ALLEEG.chanlocs);
        chanlocsName = ALLEEG(1).chaninfo.filename; % Added by Romain 18 Aug. 2010
        n_tmp   = eval(  result{1} );
        index    = eval( result{2} );
        th=result{3};
        ics=eval( result{4} );
        badcomps=result{5};
        title=result{6};
        clname=result{7};
  
        [CORRMAP,STUDY,ALLEEG] = corrmap(STUDY,ALLEEG,n_tmp,index,'chanlocs',chanlocs,'th',th,'ics',ics,'title',title,'clname',clname,'badcomps',badcomps);
        % com = sprintf('pop_corrmap(STUDY,ALLEEG,%g,%g,''chanlocs'',''%s'',''th'',''%s'',''ics'',%g,''title'',''%s'',''clname'',''%s'',''badcomps'',''%s'');', n_tmp,index,chanlocs,th,ics,title,clname,badcomps);
        com = sprintf('pop_corrmap(STUDY,ALLEEG,%g, %g,''chanlocs'',''%s'',''th'',''%s'',''ics'',%g,''title'',''%s'',''clname'',''%s'',''badcomps'',''%s'');', n_tmp,index,chanlocsName,th,ics,title,clname,badcomps);
    end

else
    % % % decode input parameters

    g = finputcheck(varargin, { 'chanlocs'  'string' [] '';...
        'th'     'string'    []     'auto' ;...
        'ics'    'integer'  [1 2 3]    2 ;....
        'title'   'string'  []     '';...
        'clname'  'string'  []    '';...
        'badcomps' 'string' {'yes','no'} 'no'});

 
    if isstr(g), error(g); end;

    if ~isempty(g.chanlocs)

        chanlocs = g.chanlocs;

    else
        chanlocs='';
    end

    if ~isempty(g.th)

        th = g.th;

    else
        th='auto';
    end

    if ~isempty(g.ics)
        if g.ics>3
        fprintf('Error "number of ICs": Maximum number allowed is 3.\n');
        fprintf('see >> help corrmap. \n');
        return
    else
        ics = g.ics;
    end
    else
        ics=2;
    end

    if ~isempty(g.title)
        title=g.title;
    else
        title='';
    end

    if ~isempty(g.clname)
        clname=g.clname;
    else
        clname='';
    end

    if ~isempty(g.badcomps)
        badcomps=g.badcomps;
    else
        badcomps='no';
    end

    [CORRMAP,STUDY,ALLEEG] = corrmap(STUDY,ALLEEG,n_tmp,index,'chanlocs',chanlocs,'th',th,'ics',ics,'title','','clname','','badcomps',badcomps);
end

