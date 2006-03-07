% cls_spec() - Returns the ICA component spectra for a dataset. Updates the EEG structure 
%              in the Matlab environment and in the .set file as well. Saves the spectra 
%              in a float file.
% Usage:    
%           >> [EEG_etc, X, f, overwrite] = cls_spec(EEG, components, ...
%                                                   freqrange, specargs, overwrite);
%
%              Computes the mean spectra of the activites of specified components of the 
%              supplied dataset. The spectra are saved in a float file. Also saves a pointer 
%              to this file in the EEG structure. If such a float file already exists, 
%              loads the spectral information from this file.  
%              Options (below) specify which components to use, and the desired frequency 
%              range. There is also an option to specify other spectopo() input variables 
%              (see >> help spectopo for details).
%
%              Returns the removed mean spectra of the selected ICA components in the 
%              requested frequency range. If the spectra were computed previously but a
%              different frequency range is selected, there is an overwrite option. 
%              so. The function will load previously computed log spectra, if any, and 
%              will remove the mean from the requested frequency range. The frequencies 
%              vector is also returned and EEG.etc is modified with the pointer to the 
%              spectra float file and with relevant information about it. 
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
%   EEG_etc   - the EEG dataset .etc structure (EEG.etc), modified with the filename and other
%               information about the float file that holds the spectra. If the spectra file 
%               already exists and wasn't modified, this output will be empty. Added/modified 
%               fields:  EEG.etc.icaspec, EEG.etc.icaspecparams, EEG.etc.icaspecmparams
%
%   X         - the mean spectra (in dB) of the requested ICA components in the selected 
%               frequency range (with the mean of each spectrum removed). 
%   f         - a vector of frequencies at which the spectra have been computed. 
%   overwrite - same as the input option, possibly modified by the user decision
%               from the pop-up menu
%
% Files output or overwritten: [dataset_filename].icaspec, 
%                              [dataset_filename].icaspecm
% 
%  See also  spectopo(), cls_erp(), cls_ersp(), cls_scalp(), eeg_preclust(), eeg_createdata()
%
% Authors:  Hilit Serby & Arnaud Delorme, SCCN, INC, UCSD, January, 2005

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

function [EEG_etc, X, f, overwrite] = cls_spec(EEG, comp, freqrange, arg ,overwrite)
    
if nargin < 1
    help cls_spec;
    return;
end;

if isfield(EEG,'icaweights')
   numc = size(EEG.icaweights,1);
else
   error('EEG.icaweights not found');
end
if ~exist('comp') | isempty(comp)
   comps = 1:numc;
else
   comps = comp; % use specified components
end
    
EEG_etc = [];
if ~exist('overwrite')
    overwrite = 0;   % default
end

if isfield(EEG,'etc')
     % if spectrum information found in datasets
     if isfield(EEG.etc, 'icaspec') & exist(fullfile(EEG.filepath, [ EEG.etc.icaspec 'm'])) & ... 
             exist(fullfile(EEG.filepath, [ EEG.etc.icaspec ]))
         params = EEG.etc.icaspecparams;
         if iscell(params)
             d= params{1};
         else
             d = params(1);
         end

         if overwrite ~= 1 % == 2: don't overwrite, read spectra info from file if exists
             mparams = EEG.etc.icaspecmparams;
             if iscell(mparams)
                 md= mparams{1};
             else
                 md = mparams(1);
             end
             fave = floatread( fullfile(EEG.filepath, [ EEG.etc.icaspec 'm']), [md 1], [], 0);
             f    = floatread( fullfile(EEG.filepath, EEG.etc.icaspec), [d 1], [], 0);

             % check whether requested information already exists
             if ~isempty('freqrange')
                 maxind = max(find(f <= freqrange(end)));
                 minind = min(find(f >= freqrange(1)));
             else
                 minind = 1;
                 maxind = md;
             end
             if ((f(maxind) == fave(md)) & (f(minind) == fave(1))) | overwrite == 2
                 X = zeros(length(comps),maxind-minind+1) ;
                 for k = 1:length(comps)
                     X(k,:) = floatread(fullfile(EEG.filepath, ...
                                        [ EEG.etc.icaspec 'm']), [md 1],[],md*comps(k))';
                 end
                 f = fave;
                 return
             end
             
             disp('Re-using existing spectrum but with new frequency boundaries');
             disp('To recompute the spectra, first delete files in the directory');
             disp('   of this dataset with extensions .icaspec and .icaspecm');

             [EEG_etc, X, f, overwrite] = cls_spec(EEG, comp, freqrange, arg, 1); % overwrite !?
             return

         else % overwrite == 1 --> overwrite existing spectra using existing spectra
              %                    but with new frequency boundaries. 

             f = floatread(fullfile(EEG.filepath, EEG.etc.icaspec), [d 1]);
             if ~isempty(freqrange)
                 maxind = max(find(f <= freqrange(end)));
                 minind = min(find(f >= freqrange(1)));
                 if f(end) < freqrange(end)
                     disp(...
['Warning! Requested high frequency limit, ' num2str(freqrange(end)) 'Hz, is out of bounds.']);
                 end
                 if f(1) > freqrange(1)
                     disp(...
['Warning! Requested low frequency limit, ' num2str(freqrange(1)) 'Hz, is out of bounds.']);
                 end
                 f = f(minind:maxind);
             else
                 minind = 1;
                 maxind = d;
             end
             X = zeros(length(comps),maxind-minind+1) ;
             for k = 1:length(comps)
                 tmp = floatread(fullfile(EEG.filepath,EEG.etc.icaspec), [d 1],[],d*comps(k));
                 X(k,:) =  tmp(minind:maxind)';
             end

             % remove the mean from each frequency across all components
             X = X - mean(X,2)*ones(1,length(f)); %remove mean
             X = X - ones(size(X,1),1)*mean(X); 

             if minind ~= 1 | maxind ~= d % new removed mean values
                 if ~isempty(arg)
                    EEG.etc.icaspecmparams = {length(f), arg{:}};
				else
                    EEG.etc.icaspecmparams = {length(f)};
				end
                %Save the updated dataset
				try
                    EEG.saved = 'no';
                    EEG = pop_saveset( EEG, 'savemode','resave');
				catch,
                    error([ 'cls_spec(): problems saving into path ' EEG.filepath])
				end
                EEG_etc = EEG.etc;

                % save removed mean spectra info in float file
                floatwrite([f X'], fullfile(EEG.filepath, ...
                                   [ EEG.filename(1:end-3) 'icaspecm']));
             end
             return
         end
     end
 end
 
% no spectra available - recompute
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
    [X, f] = spectopo(TMP.icaact, EEG.pnts, EEG.srate, 'plot', 'off', arg);  
else
    [X, f] = spectopo(TMP.icaact, EEG.pnts, EEG.srate, 'plot', 'off');
end

% save spectrum in float file
floatwrite([f X'], fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icaspec'])); 

% update the info in the dataset
EEG.etc.icaspec = [ EEG.filename(1:end-3) 'icaspec'];
if ~isempty(arg)
    EEG.etc.icaspecparams = {length(f), arg{:}};
else
    EEG.etc.icaspecparams = {length(f)};
end

% Select desired components
X = X(comps,:); 
if ~isempty('freqrange')
    maxind = max(find(f <= freqrange(end)));
    minind = min(find(f >= freqrange(1)));
    if f(end) < freqrange(end)
        disp(...
['Warning! Requested high frequency limit, ' num2str(freqrange(end)) 'Hz, is out of bounds.']);
    end
    if f(1) > freqrange(1)
        disp(...
['Warning! Requested low frequency limit, ' num2str(freqrange(1)) 'Hz, is out of bounds.']);
    end
    f = f(minind:maxind);
    X = X(:,minind:maxind);
end               

X = X - mean(X,2)*ones(1,length(f)); 
X = X - ones(size(X,1),1)*mean(X); % remove the mean from each frequency across all components
if ~isempty(arg)
    EEG.etc.icaspecmparams = {length(f), arg{:}};
else
    EEG.etc.icaspecmparams = {length(f)};
end

% save removed mean spectra info in float file
floatwrite([f X'], fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icaspecm']));

% save the updated dataset
try
    EEG.saved = 'no';
    EEG = pop_saveset( EEG, 'savemode', 'resave');
catch,
    error([ 'cls_spec: problems saving into path ' EEG.filepath])
end
EEG_etc = EEG.etc;
