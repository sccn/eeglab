% cls_scalp() - This function uses topoplot function to get the interpolated 
% Cartesian grid of the component scalp maps. The scalp map grids are saved
% into a float file and a pointer to the file is stored in the EEG structure.
% If such a float file already exists, the function loads the information from it. 
% There is an option to specify certain components. 
% The function returns the scalp map grids of all the requested components. 
% It also returns the EEG sub-structure etc (i.e EEG.etc), which is modified 
% with the pointer to the floating file and some information about the file. 
%
% Usage:
%   >> [EEG_etc, X] = cls_scalp(EEG,components);  
%   Returns the ICA scalp map grid for a dataset. 
%   Updates the EEG structure in the Matlab environment and on the disk
%   too!
%
% Inputs:
%   EEG     - an EEG data structure. 
%   components - [numeric vector] of the EEG structure for which a grid of their 
%                      scalp maps will be returned. 
%
% Outputs:
%   EEG_etc    - the EEG dataset etc structure (i.e. EEG.etc), which is
%                      modified with the pointer and some information about
%                      the floating file that holds the dataset ICA scalp map information.
%                      If the scalp map file already exists (this output will be empty). 
%   X              - the scalp map grid of the requested ICA components, each grid is 
%                     fitted into one row of X. 
%
% Authors:  Hilit Serby, Arnaud Delorme, SCCN, INC, UCSD, January, 2005
%
%  See also  topoplot, cls_scalpL, cls_scalpG, cls_erp, cls_ersp, cls_spec, eeg_preclust, eeg_createdata           

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
% Revision 1.2  2006/03/03 00:41:38  arno
% now correctly saving data
%

function [EEG_etc, X] = cls_scalp(EEG, comp)

EEG_etc = [];
%figure; toporeplot(grid,'style', 'both','plotrad', 0.5, 'intrad', 0.5, 'xsurface' ,Xi, 'ysurface',Yi );
%scalp information found in datasets
if isfield(EEG,'etc')
     if isfield(EEG.etc, 'icascalp')
         d = EEG.etc.icascalpparams; %the grid dimension 
         if iscell(d)
             d = d{1};
         end
         olddir = pwd;
         eval ( ['cd '  EEG.filepath]);
         for k = 1:length(comp)
             tmp = floatread(EEG.etc.icascalp, [d d],[],d*(d+2)*(comp(k)-1));
             tmp = tmp(:)';
             tmp = tmp(find(~isnan(tmp)));
             if k == 1
                 X = zeros(length(comp),length(tmp)) ;
             end
             X(k,:) =  tmp;
         end
         eval ([ 'cd ' olddir]); 
         return
     end
 end
 
numc = size(EEG.icaweights,1); %number of ICA comp
for k = 1:numc
    %compute scalp map grid (topoimage)
    [hfig grid plotrad Xi Yi] = topoplot( EEG.icawinv(:,k), EEG.chanlocs, 'verbose', 'off',...
        'electrodes', 'on' ,'style','both','plotrad',0.5,'intrad',0.5,'noplot', 'on', 'chaninfo', EEG.chaninfo);
    if k == 1
        d = length(grid);
        all_topos = zeros(d,(d+2)*numc);
    end
    all_topos(:,1+(k-1)*(d+2):k*(d+2) ) = [grid Yi(:,1)  Xi(1,:)'];
    %[Xi2,Yi2] = meshgrid(Yi(:,1),Xi(1,:));
end

%save topos in file
tmpfile = fullfile( EEG.filepath, [ EEG.filename(1:end-3) 'icascalp' ]); 
floatwrite(all_topos, tmpfile);
EEG.etc.icascalp       = [ EEG.filename(1:end-3) 'icascalp' ];
EEG.etc.icascalpparams = d;
try
    EEG = pop_saveset( EEG, 'savemode', 'resave');
catch,
    error([ 'cls_scalp: problems saving into path ' EEG.filepath])
end
EEG_etc = EEG.etc;

for k = 1:length(comp)
    tmp = all_topos(:,1+(k-1)*(d+2):k*(d+2)-2);
    tmp = tmp(:)';
    tmp = tmp(find(~isnan(tmp)));
    if k == 1
        X = zeros(length(comp),length(tmp)) ;
    end
    X(k,:) =  tmp;
end
