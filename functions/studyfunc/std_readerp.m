%
% std_readerp() - returns the ERP for an ICA component in an epoched dataset.
%                 The ERPs of the dataset ICA components are assumed to have 
%                 been saved in a Matlab file, [dataset_name].icaerp, in the
%                 same directory as the dataset file. If this file doesn't exist, 
%                 use std_erp() to create it, else use a pre-clustering function 
%                 that calls it: pop_preclust() or std_preclust()  
% Usage:    
%   >> [erp, times] = std_readerp(ALLEEG, setindx, component, timewindow);  
%
% Inputs:
%   ALLEEG     - an EEGLAB dataset vector (else one EEG dataset). 
%                ALLEEG must contain the dataset of interest (see 'setindx').
%   setindx    - [integer] index of the EEG dataset in the ALLEEG structure 
%                for which to read the component ERP.
%   component  - [integer] index of the component in the selected EEG dataset 
%                for which to return the ERP. 
%   timewindow - [min max] ERP time (latency) window, in ms. Must be in
%                the dataset epoch latency range.
% Outputs:
%   erp        - ERP for the requested ICA component in the selected dataset; 
%                the average of the ICA activations in all the dataset epochs.
%   times      - vector of ERP time points (latencies) in ms.
%
%  See also  std_erp(), pop_preclust(), std_preclust()          
%
% Authors: Arnaud Delorme, Hilit Serby, SCCN, INC, UCSD, February, 2005

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
% Revision 1.14  2006/03/14 01:51:02  scott
% help msg
%
% Revision 1.13  2006/03/13 23:21:27  arno
% timerange
%
% Revision 1.12  2006/03/11 07:30:41  arno
% time input
%
% Revision 1.11  2006/03/11 07:22:22  arno
% header
%
% Revision 1.10  2006/03/11 07:18:43  arno
% header information
%
% Revision 1.9  2006/03/10 16:08:19  arno
% select time range
%
% Revision 1.8  2006/03/10 00:36:34  arno
% error msg
%
% Revision 1.7  2006/03/09 18:53:15  arno
% reading all ERPs if necessary
%
% Revision 1.6  2006/03/09 18:45:40  arno
% reading all ERP
%
% Revision 1.5  2006/03/09 18:24:36  arno
% load Matlab file now
%
% Revision 1.4  2006/03/08 20:31:25  arno
% rename func
%
% Revision 1.3  2006/03/07 22:21:26  arno
% use fullfile
%
% Revision 1.2  2006/03/07 22:09:25  arno
% fix error message
%

function [X, t] = std_readerp(ALLEEG, abset, comp, timerange)

if nargin < 4
    timerange = [];
end;
    
X = [];
filename  = fullfile( ALLEEG(abset).filepath,[ ALLEEG(abset).filename(1:end-3) 'icaerp']);

for k=1:length(comp)
    try,
        erpstruct = load( '-mat', filename, [ 'comp' int2str(comp(k)) ], 'times' );
    catch
        error( [ 'Cannot read file ''' filename '''' ]);
    end;

    tmpdat    = getfield(erpstruct, [ 'comp' int2str(comp(k)) ]);
    if k == 1
        X = zeros(length(comp), length(tmpdat));
    end;
    X(k,:)  = tmpdat;
    t         = getfield(erpstruct, 'times');
end;

% select time range of interest
% -----------------------------
if ~isempty(timerange)
    maxind = max(find(t <= timerange(end)));
    minind = min(find(t >= timerange(1)));
else
    maxind = length(t);
    minind = 1;
end;
X = X(:,minind:maxind);
t = t(minind:maxind)';
