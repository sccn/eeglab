% std_spec() - Returns the ICA component spectra for a dataset. Updates the EEG structure 
%              in the Matlab environment and in the .set file as well. Saves the spectra 
%              in a file.
% Usage:    
%           >> [spec freqs] = std_spec(EEG, 'key', 'val', ...);
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
%   EEG - a loaded epoched EEG dataset structure. 
%
% Optional inputs:
%   'components' - [numeric vector] components of the EEG structure for which 
%                  activation ERPs will be computed. Note that because 
%                  computation of component spectra is relatively fast, all 
%                  components spectra are computed and saved. Only selected 
%                  component are returned by the function to Matlab
%                  {default|[] -> all}
%   'channels'   - [cell array] channels of the EEG structure for which 
%                  activation spectrum will be computed. Note that because 
%                  computation of spectrum is relatively fast, all channels 
%                  spectrum are computed and saved. Only selected channels 
%                  are returned by the function to Matlab
%                  {default|[] -> none}
%   'freqrange'  - [minhz maxhz] frequency range (in Hz) within which to 
%                  return the spectrum {default|[]: [0 sample rate/2]}. 
%   'recompute'  - ['on'|'off'] force recomputing ERP file even if it is 
%                  already on disk.
%
% Other optional spectral parameters:
%   All optional parameters to the spectopo function may be provided to this function
%   as well.
%
% Outputs:
%   spec      - the mean spectra (in dB) of the requested ICA components in the selected 
%               frequency range (with the mean of each spectrum removed). 
%   freqs     - a vector of frequencies at which the spectra have been computed. 
%
% Files output or overwritten for ICA: 
%               [dataset_filename].icaspec,   % raw spectrum of ICA components
%               [dataset_filename].icaspecm   % spectrum with the mean baseline removed
% Files output or overwritten for data: 
%               [dataset_filename].datspec, 
%               [dataset_filename].datspecm
% 
% See also  spectopo(), std_erp(), std_ersp(), std_map(), std_preclust()
%
% Authors: Arnaud Delorme, SCCN, INC, UCSD, January, 2005

% Defunct:      0 -> if frequency range is different from saved spectra, ask via a 
%                    pop-up window whether to keep existing spectra or to overwrite them. 

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
% Revision 1.33  2007/08/13 18:58:32  arno
% nothing
%
% Revision 1.32  2007/08/13 01:18:45  arno
% update help message
%
% Revision 1.31  2007/02/28 12:04:59  arno
% force recompute
%
% Revision 1.29  2006/11/21 21:59:45  arno
% computing ICA activity bug
%
% Revision 1.28  2006/11/10 00:11:12  arno
% default specmode is using spectopo
%
% Revision 1.27  2006/10/02 11:42:11  arno
% plotting scalp maps
%
% Revision 1.25  2006/03/10 22:40:21  arno
% saving average spectrum
%
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

function [X, f, overwrt] = std_spec(EEG, varargin)

overwrt = 1; % deprecated
if nargin < 1
    help std_spec;
    return;
end;

% decode inputs
% -------------
if ~isempty(varargin) 
    if ~isstr(varargin{1})
        varargin = { varargin{:} [] [] };
        if all(varargin{1} > 0) 
            options = { 'components' varargin{1} 'freqrange' varargin{2} };
        else
            options = { 'channels' -varargin{1} 'freqrange' varargin{2} };
        end;
    else
        options = varargin;
    end;
else
    options = varargin;
end;

[g spec_opt] = finputcheck(options, { 'components' 'integer' []         [];
                                      'channels'   'cell'    {}         {};
                                      'timerange'  'float'   []         [];
                                      'specmode'   'string'  {'fft' 'psd'} 'psd';
                                      'recompute'  'string'  { 'on' 'off' } 'off';
                                      'nfft'       'integer' []         [];
                                      'freqrange'  'real'    []         [] }, 'std_spec', 'ignore');
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
    filename = fullfile( EEG.filepath,[ EEG.filename(1:end-3) 'datspec']);
    prefix = 'chan';
else    
    filename = fullfile( EEG.filepath,[ EEG.filename(1:end-3) 'icaspec']);
    prefix = 'comp';
end;

% SPEC information found in datasets
% ---------------------------------
if exist(filename) & strcmpi(g.recompute, 'off')

    if strcmpi(prefix, 'comp')
        [X, f] = std_readspec(EEG, 1, g.components, g.freqrange);
    else
        [X, f] = std_readspec(EEG, 1, g.channels, g.freqrange);
    end;
    return;
    
end 
 
% No SPEC information found
% ------------------------
if isstr(EEG.data)
    TMP = eeg_checkset( EEG, 'loaddata' ); % load EEG.data and EEG.icaact
else
    TMP = EEG;
end
if strcmpi(prefix, 'comp') & isempty(TMP.icaact)
    TMP.icaact = (TMP.icaweights*TMP.icasphere)* ...
        reshape(TMP.data(TMP.icachansind,:,:), [ length(TMP.icachansind) size(TMP.data,2)*size(TMP.data,3) ]);
end;
if strcmpi(prefix, 'comp'), X = TMP.icaact;
else                        X = TMP.data;
end;

if ~isempty(g.timerange)
    timebef = find(EEG.times > g.timerange(1) & EEG.times < g.timerange(2) );
    X       = X(:,timebef,:);
end;
if strcmpi(g.specmode, 'psd')
     [X, f] = spectopo(X, EEG.pnts, EEG.srate, 'plot', 'off', 'nfft', g.nfft, spec_opt{:});  
else
    tmp   = fft(X, g.nfft, 2);
    f     = linspace(0, EEG.srate/2, size(tmp,2)/2);
    f     = f(2:end); % remove DC (match the output of PSD)
    tmp   = tmp(:,2:size(tmp,2)/2,:);
    X     = 10*log10(mean(abs(tmp).^2,3));     
end;

% Save SPECs in file (all components or channels)
% ----------------------------------
options = { spec_opt{:} 'timerange' g.timerange 'nfft' g.nfft 'specmode' g.specmode };
if strcmpi(prefix, 'comp')
    savetofile( filename, f, X, 'comp', 1:size(X,1), options);
    [X f] = std_readspec(EEG, 1, g.components, g.freqrange);
else
    savetofile( filename, f, X, 'chan', 1:size(X,1), options, { TMP.chanlocs.labels });
    [X f] = std_readspec(EEG, 1, g.channels, g.freqrange);
end;
return;

% -------------------------------------
% saving SPEC information to Matlab file
% -------------------------------------
function savetofile(filename, f, X, prefix, comps, params, labels);
    
    disp([ 'Saving SPECTRAL file ''' filename '''' ]);
    allspec.freqs      = f;
    allspec.parameters = params;
    allspec.datatype   = 'SPECTRUM';
    for k = 1:length(comps)
        allspec = setfield( allspec, [ prefix int2str(comps(k)) ], X(k,:));
    end;
    allspec.average_spec = mean(X,1);
    std_savedat(filename, allspec);

