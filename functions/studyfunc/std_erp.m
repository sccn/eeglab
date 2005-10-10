% Usage:    
%   >> [EEG_etc, data, time] = cls_erp(EEG,components,timewindow);  
%   Returns the ICA activation ERP information for a dataset. 
%   Updates the EEG structure in the Matlab environment and on the disk
%   too!
%
% cls_erp() - This function constructs the ERP of the ICA activation of a dataset,
% saves it into a float file and saves a pointer to it in the EEG structure.
% If such a float file already exists, the function loads the information. 
% There is an option to specify only certain components, and a specific time range within
% the epoch limits. The function returns the ERP of the selected ICA components in 
% the requested time window. The time window samples vectors is returned
% too, as well as the EEG sub-structure etc (i.e EEG.etc), which is modified 
% with the pointer to the floating file and some information about it. 
%
%
% Inputs:
%   EEG     - an EEG data structure. 
%   components - [numeric vector] of the EEG structure for which an ERP of their 
%                      activation will be computed. 
%   timewindow - [minms maxms] the time window limits to compute the ERP.
%
% Outputs:
%   EEG_etc    - the EEG dataset etc structure (i.e. EEG.etc), which is
%                      modified with the pointer and some information about
%                      the floating file that holds the dataset ICA ERP information.
%                      If the ERP file already exists (this output will be empty). 
%   data         - the ERP of the requested ICA components in the selected 
%                     time window. 
%   time         - a time vector of the time points in which the ERPs were computed. 
%                     The time vector is the epoch time points in the time window range. 
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
                 disp(['Warning! Requested high time window limit, ' num2str(timerange(end)) 'ms, is out of bound.']);
             end
             if t(1) > timerange(1)
                 disp(['Warning! Requested low time window limit, ' num2str(timerange(1)) 'ms, is out of bound.']);
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
    error([ 'cls_erp: problems saving into path ' EEG.filepath])
end
EEG_etc = EEG.etc;
eval ([ 'cd ' olddir]); 

%Slect desired components
X = X(comp,:); 
if ~isempty('timerange')
    maxind = max(find(EEG.times <= timerange(end)));
    minind = min(find(EEG.times >= timerange(1)));
    if EEG.times(end) < timerange(end)
         disp(['Warning! Requested high time window limit, ' num2str(timerange(end)) 'ms, is out of bound.']);
     end
     if EEG.times(1) > timerange(1)
         disp(['Warning! Requested low time window limit, ' num2str(timerange(1)) 'ms, is out of bound.']);
     end
    t = EEG.times(minind:maxind);
    X = X(:,minind:maxind);
end  
