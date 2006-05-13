% std_erp() -   Constructs and returns ICA activation ERPs for a dataset. 
%               Updates the EEG structure both in the Matlab environment 
%               and on disk. Saves the ERPs into a Matlab file, 
%               [dataset_name].icaerp, in the same directory as the dataset
%               file.  If such a file already exists, loads its information. 
%               Options allow limiting ERP coomputation to specified components
%               and to a specific latency range (within epoch limits). Returns
%               the ERP of the selected ICA components in the requested 
%               time range.
% Usage:    
%            >> [erp, times] = std_erp(EEG,components,time_range);  
% Inputs:
%   EEG          - a loaded epoched EEG dataset structure. 
%   components   - [numeric vector] components of the EEG structure for which 
%                  activation ERPs will be computed. {default|[] -> all}
%   time_range   - [minms maxms] latency window limits (in ms) within which to 
%                  compute ERPs {default|[]: [EEG.minms EEGmaxms]}
% Outputs:
%   erp         - ERP for the requested ICA components in the selected 
%                 latency window. ERPs are scaled by the RMS over of the
%                 component scalp map projection over all data channels.
%   times       - vector of times (epoch latencies in ms) for the ERP
%
% File output:     [dataset_file].icaerp     % component erp file
%
% See also:    std_spec(), std_ersp(), std_topo(), std_preclust()
%
% Authors: Arnaud Delorme, Hilit Serby, SCCN, INC, UCSD, January, 2005

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
% Revision 1.23  2006/03/14 03:28:10  scott
% help msg
%
% Revision 1.22  2006/03/11 07:08:03  arno
% header
%
% Revision 1.21  2006/03/10 16:23:43  arno
% reprogram timerange
%
% Revision 1.20  2006/03/10 00:30:21  arno
% update header
%
% Revision 1.19  2006/03/09 18:52:58  arno
% saving ERP to float file
%
% Revision 1.18  2006/03/08 22:24:37  scott
% help msg  -sm
%
% Revision 1.17  2006/03/08 22:05:19  arno
% remove bebug msg
%
% Revision 1.16  2006/03/08 22:04:47  arno
% move return to the right place
%
% Revision 1.15  2006/03/08 22:01:45  arno
% remve debug message
%
% Revision 1.14  2006/03/08 22:01:10  arno
% do not recompute for plotting
%
% Revision 1.13  2006/03/08 21:52:44  arno
% typo
%
% Revision 1.12  2006/03/08 21:51:25  arno
% fix typo
%
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

function [X, t] = std_erp(EEG, comps, timerange)

if nargin < 1
    help std_erp;
    return;
end;
    
if isfield(EEG,'icaweights')
   numc = size(EEG.icaweights,1);
else
   error('EEG.icaweights not found');
end
if nargin < 2
   comps = 1:numc;
elseif isempty(comps)
   comps = 1:numc;
end

EEG_etc = [];
if ~exist('timerange')
    timerange = [];
end

% filename 
% --------
filenameerp = fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icaerp' ]);

% ERP information found in datasets
% ---------------------------------
if exist(filenameerp)

    [X, t] = std_readerp(EEG, 1, comps, timerange);
    return;
    
end 
   
% No ERP information found
% ------------------------
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
% --------------------
if EEG.trials > 1 %epoched data
    time0 = find(EEG.times < 0);
    time0 = find(EEG.times(time0) > timerange(1));
    if ~isempty(time0)
        TMP.icaact = rmbase(TMP.icaact,EEG.pnts, time0);
    else
        TMP.icaact = rmbase(TMP.icaact,EEG.pnts);
    end
else
    TMP.icaact = rmbase(TMP.icaact);
end
TMP.icaact = reshape(TMP.icaact, [ size(TMP.icaact,1) size(TMP.data,2) size(TMP.data,3) ]);
X = repmat(sqrt(mean(TMP.icawinv.^2))', [1 TMP.pnts]) .* mean(TMP.icaact,3); % calculate ERP

% Save ERPs in file (all components)
% ----------------------------------
savetofile( filenameerp, EEG.times, X, 1:size(X,1));
[X,t] = std_readerp( EEG, 1, comps, timerange);

% -------------------------------------
% saving ERP information to Matlab file
% -------------------------------------
function savetofile(filename, t, X, comps);
    
    allerp.times      = t;
    allerp.datatype   = 'ERP';
    for k = 1:length(comps)
        allerp = setfield( allerp, [ 'comp' int2str(comps(k)) ], X(k,:));
    end;
    std_savedat(filename, allerp);
