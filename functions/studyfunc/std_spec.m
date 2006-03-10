% std_spec() - Returns the ICA component spectra for a dataset. Updates the EEG structure 
%              in the Matlab environment and in the .set file as well. Saves the spectra 
%              in a file.
% Usage:    
%           >> [EEG_etc, X, f, overwrite] = std_spec(EEG, components, ...
%                                                   freqrange, specargs, overwrite);
%
%              Computes the mean spectra of the activites of specified components of the 
%              supplied dataset. The spectra are saved in a Matlab file. If such a file 
%              already exists, loads the spectral information from this file.  
%              Options (below) specify which components to use, and the desired frequency 
%              range. There is also an option to specify other spectopo() input variables 
%              (see >> help spectopo for details).
%
%              Returns the removed mean spectra of the selected ICA components in the 
%              requested frequency range. If the spectra were computed previously but a
%              different frequency range is selected, there is an overwrite option. 
%              so. The function will load previously computed log spectra, if any, and 
%              will remove the mean from the requested frequency range. The frequencies 
%              vector is also returned. 
% Inputs:
%   EEG        - an EEGLAB dataset structure. 
%   components - [numeric vector] components of the EEG structure for which the mean 
%                spectra  are to be computed {default|[] -> all}
%   freqrange  - [minHz maxHz] the frequency range in which to compute the spectra.
%                {default: }
%   specargs   - {'key1', 'val1',...} cell array of optional spectopo inputs 
%                {default empty}
%   overwrite  - [1|2] 1 -> overwrite the saved spectra for this dataset 
%                      2 -> keep the spectra {default = 2}
% Outputs:
%   X         - the mean spectra (in dB) of the requested ICA components in the selected 
%               frequency range (with the mean of each spectrum removed). 
%   f         - a vector of frequencies at which the spectra have been computed. 
%   overwrite - same as the input option, possibly modified by the user decision
%               from the pop-up menu
%
% Files output or overwritten: 
%               [dataset_filename].icaspec, 
%               [dataset_filename].icaspecm
% 
%  See also  spectopo(), std_erp(), std_ersp(), std_map(), std_preclust()
%
% Authors:  Arnaud Delorme, Hilit Serby, SCCN, INC, UCSD, January, 2005

% Defunct:      0 -> if frequency range is different from saved spectra, ask via a 
%                    pop-up window whether to keep existing spectra or to overwrite them. 

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
% Revision 1.24  2006/03/10 17:04:56  arno
% only one spectrum version
%
% Revision 1.23  2006/03/09 18:07:45  arno
% remove output argument
%
% Revision 1.22  2006/03/09 18:05:45  arno
% reprogramming function
%
% Revision 1.21  2006/03/09 00:48:39  arno
% do not resave dataset if changing limits
%
% Revision 1.20  2006/03/09 00:39:31  arno
% erase allspec first
%
% Revision 1.19  2006/03/09 00:37:35  arno
% now writing matlab file
%
% Revision 1.18  2006/03/08 20:28:13  arno
% rename func
%
% Revision 1.17  2006/03/07 22:31:01  arno
% typo in test
%
% Revision 1.16  2006/03/07 22:28:58  arno
% fix header
%
% Revision 1.15  2006/03/07 03:56:49  scott
% edited help msg; made accept specified components -sm
%
% Revision 1.14  2006/03/06 23:16:31  arno
% change field for resave
%
% Revision 1.13  2006/03/06 23:11:52  arno
% remove comments
%
% Revision 1.12  2006/03/05 15:52:25  arno
% allowing overwrite == 2
%
% Revision 1.11  2006/03/05 00:21:20  scott
% checking on overwrite options - they look ~  -sm
%
% Revision 1.10  2006/03/04 04:11:45  scott
% edited help and msgs  -sm
%
% Revision 1.9  2006/03/03 21:42:55  arno
% new message; remove GUI
%
% Revision 1.8  2006/03/03 21:37:17  arno
% new message
%
% Revision 1.7  2006/03/03 21:32:07  arno
% same
%
% Revision 1.6  2006/03/03 21:29:52  arno
% update change dir, ICA computatation, read/write
%

function [X, f, overwrite] = std_spec(EEG, comps, freqrange, arg ,overwrite)
    
if nargin < 1
    help std_spec;
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
if nargin < 5
    overwrite = 0;   % default
end

% filename 
% --------
filenamespec = fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icaspec' ]);

if exist(filenamespec)
    [X f] = std_readspec(EEG, 1, comps, freqrange);
    return;
end
 
% no spectra available - recompute
% --------------------------------
if isstr(EEG.data)
    TMP = eeg_checkset( EEG, 'loaddata' ); % load EEG.data and EEG.icaact
else
    TMP = EEG;
end
if isempty(TMP.icaact)
    TMP.icaact = (TMP.icaweights*TMP.icasphere)* ...
                 reshape(TMP.data  , [ size(TMP.data,1)   size(TMP.data,2)*size(TMP.data,3) ]);
    TMP.icaact = reshape(TMP.icaact, [ size(TMP.icaact,1) size(TMP.data,2)*size(TMP.data,3) ]);
end;
if ~isempty(arg)
    [X, f] = spectopo(TMP.icaact, EEG.pnts, EEG.srate, 'plot', 'off', arg{:});  
else
    [X, f] = spectopo(TMP.icaact, EEG.pnts, EEG.srate, 'plot', 'off');
end

% save spectrum in file
% ---------------------
savetofile( filenamespec, f, X, 1:size(X,1), arg);
[X f] = std_readspec(EEG, 1, comps, freqrange);

% ------------------------------------------
% saving spectral information to Matlab file
% ------------------------------------------
function savetofile(filename, f, X, comps, params);
    
    allspec.freqs      = f;
    allspec.parameters = params;
    allspec.datatype   = 'SPECTRUM';
    for k = 1:length(comps)
        allspec = setfield( allspec, [ 'comp' int2str(comps(k)) ], X(k,:));
    end;
    allspec.average_spec = mean(X,1);
    std_savedat(filename, allspec);
