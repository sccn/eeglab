%   cls_erp() - Returns ICA activation ERP information for a dataset. 
%               Updates the EEG structure both in the Matlab environment 
%               and on disk. Constructs the ERP of the ICA activation of a dataset,
%               saves it into a float file, and then saves a pointer to it 
%               in the EEG structure. If such a float file already exists, 
%               the function loads its information. There is an option to 
%               limit ERP coomputation to specified components, and to a 
%               specific latency range within the epoch limits. The function 
%               returns the ERP of the selected ICA components in the requested 
%               time window. A latencies vector is returned too, as well as 
%               the EEG sub-structure .etc (i.e., EEG.etc), to which is added 
%               a pointer to the ERP float file and some information about it. 
% Usage:    
%        >> [EEG_etc, data, times] = cls_erp(EEG,components,timewindow);  
%
% Inputs:
%   EEG          - a loaded EEG data structure. 
%   components   - [numeric vector] components of the EEG structure for which 
%                      an activation ERP will be computed. 
%   timewindow   - [minms maxms] the latency window limits within which to compute 
%                      the ERP.
% Outputs:
%   EEG_etc      - the EEG dataset .etc sub-structure (i.e., EEG.etc), to which 
%                      is added an ERP file pointer plus some information about
%                      the float file that holds the component ERP information.
%                      If the ERP file already exists, this output will be empty. 
%   data         - ERP for the requested ICA components in the selected 
%                     latency window. 
%   times        - a vector of epoch latencies at which the ERPs are computed. 
%
%  See also  cls_spec, cls_ersp, cls_scalp, eeg_preclust, eeg_createdata           
%
% Authors:  Hilit Serby, SCCN, INC, UCSD, January, 2005

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

function [EEG_etc, X, t] = cls_erp(EEG, comp, timerange)

EEG_etc = [];
if ~exist('timerange')
    timerange = [];
end
%erp information found in datasets
if isfield(EEG,'etc')
     if isfield(EEG.etc, 'icaerp')
         d = EEG.etc.icaerpparams;
         olddir = pwd;
         eval ( ['cd '  EEG.filepath]);
         t = floatread(EEG.etc.icaerp, [d 1]);
         blind = max(find(t <= 0)); %the baseline index
         if ~isempty(timerange)
             maxind = max(find(t <= timerange(end)));
             minind = min(find(t >= timerange(1)));
             if t(end) < timerange(end)
                 disp(['Warning! The requested max latency window limit, ' num2str(timerange(end)) 'ms, is out of bounds.']);
             end
             if t(1) > timerange(1)
                 disp(['Warning! The requested min latency window limit, ' num2str(timerange(1)) 'ms, is out of bounds.']);
             end
             t = t(minind:maxind);
         else
             minind = 1;
             maxind = d;
         end
         X = zeros(length(comp),maxind-minind+1) ;
         for k = 1:length(comp)
             tmp = floatread(EEG.etc.icaerp, [d 1],[],d*comp(k));
             X(k,:) =  tmp(minind:maxind)';
         end
         eval ([ 'cd ' olddir]); 
         return
     end
 end
 
%no erp information
if isempty(EEG.icaact)
    EEG = eeg_checkset( EEG, 'loaddata' ); %load EEG.data and EEG.icaact
end
%remove base line
if EEG.trials > 1 %epoched data
    time0 = find(EEG.times==0);
    if ~isempty(time0)
        tmp = rmbase(EEG.icaact,EEG.pnts,1:time0);
    else
        tmp = rmbase(EEG.icaact);
    end
else
    tmp = rmbase(EEG.icaact);
end
tmp = reshape(tmp,EEG.nbchan,EEG.pnts,EEG.trials);
X = mean(tmp,3); %calculate ERP
t = EEG.times';
%save erp in file
olddir = pwd;
eval ( ['cd '  EEG.filepath]);
floatwrite([t X'], [ EEG.filename(1:end-3) 'icaerp']);
%update the info in the dataset
EEG.etc.icaerp = [ EEG.filename(1:end-3) 'icaerp'];
EEG.etc.icaerpparams = length(t);
try
    EEG = pop_saveset( EEG, 'filename', EEG.filename, 'filepath', EEG.filepath, 'savemode','twofiles');
catch,
    error([ 'problem saving information into path ' EEG.filepath])
end
EEG_etc = EEG.etc;
eval ([ 'cd ' olddir]); 

%Slect desired components
X = X(comp,:); 
if ~isempty('timerange')
    maxind = max(find(EEG.times <= timerange(end)));
    minind = min(find(EEG.times >= timerange(1)));
    if EEG.times(end) < timerange(end)
         disp(['Warning! Requested max latency window limit, ' num2str(timerange(end)) 'ms, is out of bounds.']);
     end
     if EEG.times(1) > timerange(1)
         disp(['Warning! Requested min latency window limit, ' num2str(timerange(1)) 'ms, is out of bounds.']);
     end
    t = EEG.times(minind:maxind);
    X = X(:,minind:maxind);
end  
