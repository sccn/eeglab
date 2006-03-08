%   std_erp() - Constructs and returns ICA activation ERPs for a dataset. 
%               Updates the EEG structure both in the Matlab environment 
%               and on disk Saves the ERPs into a float file, and then saves 
%               a pointer to the file in the EEG structure. If such a float file 
%               already exists, loads its information. Options allow
%               limiting ERP coomputation to specified components, and to a 
%               specific latency range within the epoch limits. The function 
%               returns the ERP of the selected ICA components in the requested 
%               time window. A latencies vector is returned, as well as the
%               EEG sub-structure EEG.etc, to which is added a pointer to the 
%               ERP float file and some information about it. 
% Usage:    
%            >> [EEG_etc, data, times] = std_erp(EEG,components,timewindow);  
% Inputs:
%   EEG          - a loaded epoched EEG dataset structure. 
%   components   - [numeric vector] components of the EEG structure for which 
%                      activation ERPs will be computed. {default|[] -> all}
%   timewindow   - [minms maxms] latency window limits (in ms) within which to 
%                      compute ERPs {default|[]: [EEG.minms EEGmaxms]}
% Outputs:
%   EEG_etc      - the EEG dataset .etc sub-structure (i.e., EEG.etc), to which 
%                      is added an ERP file pointer plus some information about
%                      the float file that holds the component ERP information.
%                      If the ERP file already exists, this output will be empty. 
%   data         - ERP for the requested ICA components in the selected 
%                     latency window. 
%   times        - a vector of epoch latencies at which the ERPs are computed. 
%
% File output:     [dataset_file].icaerp
%                  [dataset_file].set
%
% See also  std_spec(), std_ersp(), std_topo(), std_preclust()
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

% $Log: not supported by cvs2svn $
% Revision 1.11  2006/03/08 20:29:26  arno
% rename func
%
% Revision 1.10  2006/03/07 22:40:10  arno
% floatwrite in double
%
% Revision 1.9  2006/03/07 03:27:25  scott
% accepting [] component list -sm
%
% Revision 1.8  2006/03/07 03:24:05  scott
% reworked help msg; clarified filename output; made function accept default comps
% -sm
%
% Revision 1.7  2006/03/06 23:17:09  arno
% change fields for resave
%
% Revision 1.6  2006/03/03 23:34:18  arno
% recomputing if our of bound
%
% Revision 1.5  2006/03/03 22:58:55  arno
% update call to pop_saveset
%
% Revision 1.4  2006/03/03 22:51:25  arno
% fix same thing
%
% Revision 1.3  2006/03/03 22:44:54  arno
% [6~[6~floatread/flotwrite folder fix; computation of ICA fix
%

function [EEG_etc, X, t] = std_erp(EEG, comp, timerange)

if nargin < 1
    help std_erp;
    return;
end;
    
if ~exist('comp') | isempty(comp)  
   comps = 1:size(EEG.icaweights,1);
else
   comps = comp;
end

EEG_etc = [];
if ~exist('timerange')
    timerange = [];
end

% ERP information found in datasets
if isfield(EEG,'etc')
     if isfield(EEG.etc, 'icaerp') & exist(fullfile(EEG.filepath, EEG.etc.icaerp))
         d = EEG.etc.icaerpparams;
         t = floatread(fullfile(EEG.filepath, EEG.etc.icaerp), [d 1]);
         blind = max(find(t <= 0)); %the baseline index
         if ~isempty(timerange)
             maxind = max(find(t <= timerange(end)));
             minind = min(find(t >= timerange(1)));
             if t(end) < timerange(end) | t(1) > timerange(1)
                 disp(['Warning! The requested max latency window limit is outside the previous bounds; recomputing the ERPs...']);
                 EEG.etc = rmfield(EEG.etc, 'icaerp');
                 [EEG_etc, X, t] = std_erp(EEG, comp, timerange);
                 return;
             end
             t = t(minind:maxind);
         else
             minind = 1;
             maxind = d;
         end
         X = zeros(length(comps),maxind-minind+1) ;
         for k = 1:length(comps)
             tmp = floatread(fullfile(EEG.filepath, EEG.etc.icaerp), [d 1],[],d*comps(k));
             X(k,:) =  tmp(minind:maxind)';
         end

         %save erp in file
         % ---------------
         floatwrite(double([t X']), fullfile( EEG.filepath, [ EEG.filename(1:end-3) 'icaerp']));
         
         %update the info in the dataset
         % -----------------------------
         EEG.etc.icaerp       = [ EEG.filename(1:end-3) 'icaerp'];
         EEG.etc.icaerpparams = length(t);
         try
             EEg.saved = 'no';
             EEG = pop_saveset( EEG, 'savemode','resave');
         catch,
             error([ 'problem saving information into path ' EEG.filepath])
         end
         EEG_etc = EEG.etc;
         return;

     end
 end
  
% No ERP information found
if isstr(EEG.data)
    TMP = eeg_checkset( EEG, 'loaddata' ); % load EEG.data and EEG.icaact
else
    TMP = EEG;
end
if isempty(TMP.icaact)
    TMP.icaact = (TMP.icaweights*TMP.icasphere)* ...
                 reshape(TMP.data  , [ size(TMP.data,1)   size(TMP.data,2)*size(TMP.data,3) ]);
    TMP.icaact = reshape(TMP.icaact, [ size(TMP.icaact,1) size(TMP.data,2) size(TMP.data,3) ]);
end;

% Remove baseline mean
if EEG.trials > 1 %epoched data
    time0 = find(EEG.times < 0);
    if ~isempty(time0)
        TMP.icaact = rmbase(TMP.icaact,EEG.pnts, time0);
    else
        TMP.icaact = rmbase(TMP.icaact,EEG.pnts);
    end
else
    TMP.icaact = rmbase(TMP.icaact);
end
maxind = max(find(EEG.times <= timerange(end)));
minind = min(find(EEG.times >= timerange(1)));
TMP.icaact = reshape(TMP.icaact, [ size(TMP.icaact,1) size(TMP.data,2) size(TMP.data,3) ]);
X = mean(TMP.icaact(:,minind:maxind,:),3); % calculate ERP
t = EEG.times(minind:maxind)';

% Save ERPs in file
% ---------------
floatwrite(double([t X']), fullfile( EEG.filepath, [ EEG.filename(1:end-3) 'icaerp']));

% Update the ERP info in the dataset .etc sub-structure
% -----------------------------------------------------
EEG.etc.icaerp       = [ EEG.filename(1:end-3) 'icaerp'];
EEG.etc.icaerpparams = length(t);
try
    EEG.saved = 'no';
    EEG = pop_saveset( EEG, 'savemode','resave');
catch,
    error([ 'problem saving information into path ' EEG.filepath])
end
EEG_etc = EEG.etc;

% Select desired components
X = X(comps,:); 
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
