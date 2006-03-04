%
% cls_spec() - Returns the ICA component spectra for a dataset. Updates the EEG structure in the 
%              Matlab environment and in the .set file as well. Saves the spectra in
% Usage:    
%   >> [EEG_etc, X, f, overwrite] = cls_spec(EEG, components, freqrange, specargs, overwrite)
%
%              This function computes the mean spectra of the activites of specfied 
%              components of the supplied dataset. The spectra are saved in a float file. 
%              Also saves a pointer to this file in the EEG structure. If such a float file 
%              already exists, the function loads the spectral information from this file.  
%              Options (below) specify which components to use, and the desired frequency 
%              range. There is also an option to set different spectopo input variables 
%              (see >> help spectopo for details).
%
%              cls_spec() returns the removed mean spectra of the selected ICA components 
%              in the requested frequency range. If the spectra were computed previously
%              but a different frequency range is selected, there is an overwrite variable 
%              to do so. Depending on the variable setting, this can open a pop-up query 
%              window to ask the user for permission to overwrite the previous spectra.
%              Otherwise, the function will load the pre-computed spectra, and will remove 
%              the mean from the requested frequency range. The frequency vector is returned
%              as well, and the EEG sub-structure EEG.etc is modified with the pointer to 
%              the spectra float file and with some information about it. 
% Inputs:
%   EEG        - an EEGLAB dataset structure. 
%   components - [numeric vector] components of the EEG structure for which the mean 
%                spectra  are to be computed. 
%   freqrange  - [minHz maxHz] the frequency range in which to compute the spectra.
%   specargs   - {'key1', 'val1',...} cell array of optional spectopo inputs (default empty).
%
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
% Files output or overwritten: EEG.filename with extension changed to .icaspec, .icaspecm
% 
%  See also  spectopo(), cls_erp(), cls_ersp(), cls_scalp(), eeg_preclust(), eeg_createdata()
%
% Authors:  Hilit Serby, SCCN, INC, UCSD, January, 2005

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

function [EEG_etc, X, f, overwrite] = cls_spec(EEG, comp, freqrange, arg,overwrite)
    
if nargin < 1
    help cls_spec;
    return;
end;
    
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

         if overwrite ~= 1 % = 0 or 2: don't overwrite, read spectra info from file if exists
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
             if (f(maxind) == fave(md)) & (f(minind) == fave(1))
                 X = zeros(length(comp),maxind-minind+1) ;
                 for k = 1:length(comp)
                     X(k,:) = floatread(fullfile(EEG.filepath, ...
                                        [ EEG.etc.icaspec 'm']), [md 1],[],md*comp(k))';
                 end
                 f = fave;
                 return
             end

% $$$              
% $$$              if overwrite ~= 2
% $$$                  set_yes =  [ 'set(findobj(''parent'', gcbf, ''tag'', ''spec_yesall''), ''value'', 0);' ...
% $$$                     'set(findobj(''parent'', gcbf, ''tag'', ''spec_no''), ''value'', 0);'  ...
% $$$                     'set(findobj(''parent'', gcbf, ''tag'', ''spec_noall''), ''value'', 0);'];
% $$$                  set_yesall =  [ 'set(findobj(''parent'', gcbf, ''tag'', ''spec_yes''), ''value'', 0);' ...
% $$$                     'set(findobj(''parent'', gcbf, ''tag'', ''spec_no''), ''value'', 0);'  ...
% $$$                     'set(findobj(''parent'', gcbf, ''tag'', ''spec_noall''), ''value'', 0);'];
% $$$                  set_no =  [ 'set(findobj(''parent'', gcbf, ''tag'', ''spec_yesall''), ''value'', 0);' ...
% $$$                     'set(findobj(''parent'', gcbf, ''tag'', ''spec_yes''), ''value'', 0);' ...
% $$$                     'set(findobj(''parent'', gcbf, ''tag'', ''spec_noall''), ''value'', 0);'];
% $$$                  set_noall =  [ 'set(findobj(''parent'', gcbf, ''tag'', ''spec_yesall''), ''value'', 0);' ...
% $$$                     'set(findobj(''parent'', gcbf, ''tag'', ''spec_yes''), ''value'', 0);' ...
% $$$                     'set(findobj(''parent'', gcbf, ''tag'', ''spec_no''), ''value'', 0);'  ];
% $$$                  spec_ans = inputgui({[1] [1] [1 1 ] [1 1] [1]}, ...
% $$$                          { {'style' 'text' 'string' ['Spectrum infomation [' num2str(round(fave(1))) ' ' num2str(round(fave(end)))  '] Hz already exists for dataset: '  EEG.setname '. ' ] } ...
% $$$                          {'style' 'text' 'string' 'Would you like to recalculate the spectrum and overwrite these information?' } ...
% $$$                          {'style' 'checkbox' 'tag' 'spec_yes' 'string' 'Yes' 'value' 1 'Callback' set_yes }  ...
% $$$                          {'style' 'checkbox' 'tag' 'spec_yesall' 'string' 'Yes to all datasets' 'value' 0 'Callback' set_yesall }  ...
% $$$                          {'style' 'checkbox' 'tag' 'spec_no' 'string' 'Use existing spectrum info' 'value' 0 'Callback' set_no } ...
% $$$                          {'style' 'checkbox' 'tag' 'spec_noall' 'string' 'Use existing info for all sets' 'value' 0 'Callback' set_noall } {} }, ...
% $$$                      '', 'Recalculate spectrum information -- part of cls_spec()'); 
% $$$                  switch (find(celltomat(spec_ans)))
% $$$                      case 1
% $$$                          [EEG_etc, X, f, overwrite] = cls_spec(EEG, comp, freqrange, arg, 1); % overwrite the info in this dataset
% $$$                          overwrite = 0; %but ask before overwriting the rest of the datasets
% $$$                          return;
% $$$                      case 2
% $$$                          [EEG_etc, X, f, overwrite] = cls_spec(EEG, comp, freqrange, arg, 1);
% $$$                          return;
% $$$                      case 3
% $$$                          overwrite = 0; % don't overwrite the info in this datase but keep asking for the other datasets
% $$$                      case 4
% $$$                          overwrite = 2; % keep the info in all datasets
% $$$                  end
% $$$              end
% $$$              
             
             disp('Using existing spectra but with new frequency boundaries');
             disp('To recompute spectra, first delete files in this dataset dir');
             disp(' with extensions .icaspec and .icaspecm');

             [EEG_etc, X, f, overwrite] = cls_spec(EEG, comp, freqrange, arg, 1);
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
             X = zeros(length(comp),maxind-minind+1) ;
             for k = 1:length(comp)
                 tmp = floatread(fullfile(EEG.filepath,EEG.etc.icaspec), [d 1],[],d*comp(k));
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
X = X(comp,:); 
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
    EEG = pop_saveset( EEG, 'savemode', 'resave');
catch,
    error([ 'cls_spec: problems saving into path ' EEG.filepath])
end
EEG_etc = EEG.etc;
