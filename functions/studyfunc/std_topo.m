% std_topo() - uses topoplot() to get the interpolated Cartesian grid of the 
%               specified component topo maps. The topo map grids are saved
%               into a (.icatopo) float file and a pointer to the file is stored 
%               in the EEG structure. If such a float file already exists, 
%               loads the information from it. 
%
%               Returns the topo map grids of all the requested components. Also
%               returns the EEG sub-structure etc (i.e EEG.etc), which is modified 
%               with a pointer to the float file and some information about the file. 
% Usage:
%               >> [EEG_etc, X] = std_topo(EEG, components, option);  
%
%                  % Returns the ICA topo map grid for a dataset. 
%                  % Updates the EEG structure in the Matlab environment and re-saves
% Inputs:
%   EEG        - an EEG dataset structure. 
%   components - [numeric vector] components in the EEG structure to compute topo maps
%                      {default|[] -> all}      
%   option     - ['gradient'|'laplacian'|'none'] compute gradient or laplacian of
%                the scale topography. {default is 'none' = the interpolated topo map}
% Outputs:
%   EEG_etc    - the EEG dataset EEG.etc structure modified with the file name 
%                     and information about the float file that holds the dataset 
%                     component topo maps. If the topo map file already exists 
%                     this output will be empty. 
%   X          - the topo map grid of the requested ICA components, each grid is 
%                     one ROW of X. 
%
% File output: [dataset_name].icatopo
%  
% Authors:  Hilit Serby, Arnaud Delorme, SCCN, INC, UCSD, January, 2005
%
% See also  topoplot(), std_erp(), std_ersp(), std_spec(), std_preclust()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
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
% Revision 1.7  2006/03/07 03:46:43  scott
% reworked help msg; made function accept specified component list or [] -sm
%
% Revision 1.6  2006/03/06 23:45:28  arno
% adding gradient and laplacian
%
% Revision 1.5  2006/03/06 23:17:34  arno
% fix resave
%
% Revision 1.4  2006/03/03 21:22:07  arno
% remve change folder
%
% Revision 1.3  2006/03/03 21:19:16  arno
% change topo map file name
%
% Revision 1.2  2006/03/03 00:41:38  arno
% now correctly saving data
%

function [EEG_etc, X] = std_topo(EEG, comp, option)

if nargin < 1
    help std_topo;
    return;
end;
if isfield(EEG,'icaweights')
     numc = size(EEG.icaweights,1); % number of ICA comps
else
     error('EEG.icaweights not found')
end
if nargin < 2 | isempty(comp)
    comps = 1:numc;
else
    comps = comp;
end
if nargin < 3
    option = 'none';
end;

EEG_etc = [];
% figure; toporeplot(grid,'style', 'both','plotrad', 0.5, 'intrad', 0.5, 'xsurface' ,Xi, 'ysurface',Yi );

% Topo information found in dataset
if isfield(EEG,'etc')
     if isfield(EEG.etc, 'icatopo') & exist(fullfile(EEG.filepath, EEG.etc.icatopo))
         d = EEG.etc.icatopoparams; %the grid dimension 
         if iscell(d)
             d = d{1};
         end
         for k = 1:length(comps)
             tmp = floatread(fullfile(EEG.filepath, EEG.etc.icatopo), ...
                                                        [d d],[],d*(d+2)*(comps(k)-1));
             if strcmpi(option, 'gradient')
                 [tmpx, tmpy]  = gradient(tmp); %Gradient
                 tmp = [tmpx(:); tmpy(:)]';
             elseif strcmpi(option, 'laplacian')
                 tmp = del2(tmp); %Laplacian
                 tmp = tmp(:)';
             else
                 tmp = tmp(:)';
             end;

             tmp = tmp(find(~isnan(tmp)));
             if k == 1
                 X = zeros(length(comps),length(tmp)) ;
             end
             X(k,:) =  tmp;
         end
         return
     end
 end
 
for k = 1:numc

    % compute topo map grid (topoimage)
    [hfig grid plotrad Xi Yi] = topoplot( EEG.icawinv(:,k), EEG.chanlocs, ...
                                          'verbose', 'off',...
                                           'electrodes', 'on' ,'style','both',...
                                           'plotrad',0.5,'intrad',0.5,...
                                           'noplot', 'on', 'chaninfo', EEG.chaninfo);
    if k == 1
        d = length(grid);
        all_topos = zeros(d,(d+2)*length(comps));
    end
    all_topos(:,1+(k-1)*(d+2):k*(d+2) ) = [grid Yi(:,1)  Xi(1,:)'];
    % [Xi2,Yi2] = meshgrid(Yi(:,1),Xi(1,:));
end

% Save topos in file
tmpfile = fullfile( EEG.filepath, [ EEG.filename(1:end-3) 'icatopo' ]); 
floatwrite(all_topos, tmpfile);
EEG.etc.icatopo       = [ EEG.filename(1:end-3) 'icatopo' ];
EEG.etc.icatopoparams = d;
try
    EEG.saved = 'no';
    EEG = pop_saveset( EEG, 'savemode', 'resave');
catch,
    error([ 'std_topo: problems saving into path ' EEG.filepath])
end
EEG_etc = EEG.etc;

for k = 1:length(comps)
    tmp = all_topos(:,1+(k-1)*(d+2):k*(d+2)-2);
    
    if strcmpi(option, 'gradient')
        [tmpx, tmpy]  = gradient(tmp); % Gradient
        tmp = [tmpx(:); tmpy(:)]';
    elseif strcmpi(option, 'laplacian')
        tmp = del2(tmp); % Laplacian
        tmp = tmp(:)';
    else
        tmp = tmp(:)';
    end;

    tmp = tmp(find(~isnan(tmp)));
    if k == 1
        X = zeros(length(comps),length(tmp)) ;
    end
    X(k,:) =  tmp;
end
