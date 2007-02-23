% std_erp() -   Constructs and returns ICA component activation or channel activity 
%               average ERP for dataset epochs. Saves the ERPs into a Matlab file, 
%                          '[dataset_name].icaerp'
%               for ICA components or 
%                          '[dataset_name].daterp' 
%               for scalp channels, in the same directory as the dataset file.  
%               If such a file already exists, std_erp() loads its information. 
% Usage:    
%               >> [ERP, ERPtimes] = std_erp(EEG, 'key', 'val', ...);
% Inputs:
%   EEG          - a loaded, epoched EEG dataset structure. 
%
% Optional inputs:
%   'components' - [numeric vector] components of the EEG structure for which 
%                  activation ERPs will be computed. Note that because 
%                  computation of ERPs is so fast, all components ERP are
%                  computed and saved, though only the selected components
%                  are returned by the function to the Matlab work space
%                  {default|[] -> all}
%   'channels'   - [cell array] channels of the EEG structure for which 
%                  ERPs will be computed. Note that be epochscause computation of ERPs 
%                  is so fast, all channel ERPs are computed and saved. Only 
%                  the selected channels are returned by the function to 
%                  the Matlab work space {default|[] -> all}
%   'time_range' - [minms maxms] latency window limits (in ms) within which to 
%                  compute ERPs It is not advisable to change this range unless 
%                  you require a specific baseline. The plotting functions
%                  can restrict the window plotted later.
%                  {default|[]: [EEG.min EEG.max]}. 
% Outputs:
%   ERP         - (comps/chans,times) ERP matrix for the requested ICA components 
%                 or scalp channels in the selected latency window. Component ERPs 
%                 are scaled by the RMS component scalp map projection to all the 
%                 scalp channels; they are thus in RMS microvolts. Channel ERPs are 
%                 in microvolts.
%   ERPtimes    - vector of ERP epoch latencies (in ms) 
%
% File output:     [dataset_file].icaerp     % component ERP Matlab file
%      or          [dataset_file].daterp     % channel ERP Matlab file
%                  
% Matlab variables in the file output:
%                               compN        % (or  chanN) for N = all comp or chan indices
%                               times        % epoch latencies
%                               datafile     % full datafile pathname
%                               datatype     % 'ERP'
%
% To load the component file output:   
%                               >> load -mat [dataset_file].icaerp % or std_readerp()
%
% See also:    std_erpplot(), std_readerp(), std_spec(), std_ersp(), std_topo(), std_preclust()
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, January, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, October 11, 2004, arno@sccn.ucsd.edu
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
% Revision 1.28  2007/02/20 14:21:40  arno
% option to force recomputing
%
% Revision 1.27  2007/02/20 01:46:30  scott
% completed and clarified the help msg
%
% Revision 1.26  2006/10/02 11:40:26  arno
% minor things
%
% Revision 1.24  2006/05/13 11:59:58  arno
% now clustering using RMS
%
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

function [X, t] = std_erp(EEG, varargin); %comps, timerange)

if nargin < 1
    help std_erp;
    return;
end;

% decode inputs
% -------------
if ~isempty(varargin) 
    if ~isstr(varargin{1})
        varargin = { varargin{:} [] [] };
        if all(varargin{1} > 0) 
            options = { 'components' varargin{1} 'timerange' varargin{2} };
        else
            options = { 'channels' -varargin{1} 'timerange' varargin{2} };
        end;
    else
        options = varargin;
    end;
else
    options = varargin;
end;

g = finputcheck(options, { 'components' 'integer' []         [];
                           'channels'   'cell'    {}         {};
                           'recompute'  'string'  { 'on' 'off' } 'off';
                           'timerange'  'real'    []         [EEG.xmin EEG.xmax]*1000 }, 'std_erp');
if isstr(g), error(g); end;
if isfield(EEG,'icaweights')
   numc = size(EEG.icaweights,1);
else
   error('EEG.icaweights not found');
end
if isempty(g.components)
    g.components = 1:numc;
end

EEG_etc = [];

% filename 
% --------
if ~isempty(g.channels)
    filename = fullfile( EEG.filepath,[ EEG.filename(1:end-3) 'daterp']);
    prefix = 'chan';
else    
    filename = fullfile( EEG.filepath,[ EEG.filename(1:end-3) 'icaerp']);
    prefix = 'comp';
end;

% ERP information found in datasets
% ---------------------------------
if exist(filename) & strcmpi(g.recompute, 'off')

    if strcmpi(prefix, 'comp')
        [X, t] = std_readerp(EEG, 1, g.components, g.timerange);
    else
        [X, t] = std_readerp(EEG, 1, g.channels, g.timerange);
    end;
    return;
    
end 
   
% No ERP information found
% ------------------------
if isstr(EEG.data)
    TMP = eeg_checkset( EEG, 'loaddata' ); % load EEG.data and EEG.icaact
else
    TMP = EEG;
end
if strcmpi(prefix, 'comp') & isempty(TMP.icaact)
    TMP.icaact = (TMP.icaweights*TMP.icasphere)* ...
        reshape(TMP.data(TMP.icachansind,:,:), [ length(TMP.icachansind) size(TMP.data,2)*size(TMP.data,3) ]);
    TMP.icaact = reshape(TMP.icaact, [ size(TMP.icaact,1) size(TMP.data,2) size(TMP.data,3) ]);
end;
if strcmpi(prefix, 'comp'), X = TMP.icaact;
else                        X = TMP.data;
end;

% Remove baseline mean
% --------------------
if EEG.trials > 1 %epoched data
    time0 = find(EEG.times < 0);
    time0 = find(EEG.times(time0) > g.timerange(1));
    if ~isempty(time0)
        X = rmbase(X,EEG.pnts, time0);
    else
        X = rmbase(X,EEG.pnts);
    end
else
    X = rmbase(X);
end
X = reshape(X, [ size(X,1) size(TMP.data,2) size(TMP.data,3) ]);
if strcmpi(prefix, 'comp')
    X = repmat(sqrt(mean(TMP.icawinv.^2))', [1 TMP.pnts]) .* mean(X,3); % calculate ERP
else    
    X = mean(X, 3);
end;

% Save ERPs in file (all components or channels)
% ----------------------------------
if strcmpi(prefix, 'comp')
    savetofile( filename, EEG.times, X, 'comp', 1:size(X,1));
    [X,t] = std_readerp( EEG, 1, g.components, g.timerange);
else
    savetofile( filename, EEG.times, X, 'chan', 1:size(X,1), { TMP.chanlocs.labels });
    [X,t] = std_readerp( EEG, 1, g.channels, g.timerange);
end;

% -------------------------------------
% saving ERP information to Matlab file
% -------------------------------------
function savetofile(filename, t, X, prefix, comps, labels);
    
    disp([ 'Saving ERP file ''' filename '''' ]);
    allerp.times      = t;
    allerp.datatype   = 'ERP';
    for k = 1:length(comps)
        allerp = setfield( allerp, [ prefix int2str(comps(k)) ], X(k,:));
    end;
    if nargin > 5
        allerp.labels = labels;
    end;
    std_savedat(filename, allerp);
