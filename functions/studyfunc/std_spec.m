% cls_spec() - Returns the ICA spectra of a dataset. Updates the EEG structure in the 
%              Matlab environment and on the disk too!
%
% Usage:    
%   >> [EEG_etc, X, f, overwrite] = cls_spec(EEG, components, freqrange, arg, overwrite)
%   
%
% cls_spec() - This function computes the spectra of a dataset ICA components,
% saves it into a float file and saves a pointer to it in the EEG structure.
% If such a float file already exists, the function loads the spectra information. 
% There is an option to specify only certain components, and a desired
% frequency range. There is also an option to set different spectopo input
% variables (see spectopo for details).
% The function returns the removed mean spectra of the selected ICA components in 
% the requested frequency range. If the spectrum were already computed before, 
% and a different frequency range is desired there is an overwrite variable to do so,
% which depending on the variable can open a pop-up window to ask the user. 
% This will load the pre-computed spectra, and will remove the mean from
% the requested frequency range. The frequency samples vectors is returned
% too, as well as the EEG sub-structure etc (i.e. EEG.etc), which is modified 
% with the pointer to the floating file and some information about it. 
%
%
% Inputs:
%   EEG        - an EEGLAB data structure. 
%   components - [numeric vector] of the EEG structure for which a spectrum  
%                will be computed. 
%   freqrange  - [minHz maxHz] the frequency range to compute the spectra.
%   arg        - {'key1', 'val1',...} cell array with optional spectopo inputs (default empty).
%   overwrite  - [0|1|2] 0 - if frequency range is different from saved info ask using a 
%                pop-up menu if to keep existing spectra or overwrite it, 1- overwrite
%                the saved spectra of this dataset or other datasets if exist, 2 - keep the 
%                spectra for all the datasets in STUDY (default - 0).
%
%
% Outputs:
%   EEG_etc   - the EEG dataset etc structure (i.e. EEG.etc), which is
%               modified with the pointer and some information about
%               the floating file that holds the dataset spectra information.
%               If the spectra file already exists and wasn't modified (this output will be empty). 
%   X         - the spectra of the requested ICA components in the selected 
%               frequency range (spectra mean removed). 
%   f         - a frequency vector of the points in which the spectra were computed. 
%   overwrite - same as input option, only modified with what the user
%               asked for the rest of the datasets in STUDY (from the
%               pop-up menu, or from command line).
%
%  See also  spectopo, cls_erp, cls_ersp, cls_scalp, eeg_preclust, eeg_createdata         
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

function [EEG_etc, X, f, overwrite] = cls_spec(EEG, comp, freqrange, arg,overwrite)
EEG_etc = [];
if ~exist('overwrite')
    overwrite = 0;
end
if isfield(EEG,'etc')
     if isfield(EEG.etc, 'icaspec') %spectrum information found in datasets
         params = EEG.etc.icaspecparams;
         if iscell(params)
             d= params{1};
         else
             d = params(1);
         end
         if overwrite ~= 1 % don't overwrite spectrum info if exists
             mparams = EEG.etc.icaspecmparams;
             if iscell(mparams)
                 md= mparams{1};
             else
                 md = mparams(1);
             end
             olddir = pwd;
             eval ( ['cd '  EEG.filepath]);
             fave = floatread([ EEG.etc.icaspec 'm'], [md 1], [], 0);
             f = floatread(EEG.etc.icaspec, [d 1], [], 0);
             %check if requested information already exists
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
                     X(k,:) = floatread([ EEG.etc.icaspec 'm'], [md 1],[],md*comp(k))';
                 end
                 f = fave;
                 eval ( ['cd '  olddir]);
                 return
             end
             
             if overwrite ~= 2
                 set_yes =  [ 'set(findobj(''parent'', gcbf, ''tag'', ''spec_yesall''), ''value'', 0);' ...
                    'set(findobj(''parent'', gcbf, ''tag'', ''spec_no''), ''value'', 0);'  ...
                    'set(findobj(''parent'', gcbf, ''tag'', ''spec_noall''), ''value'', 0);'];
                 set_yesall =  [ 'set(findobj(''parent'', gcbf, ''tag'', ''spec_yes''), ''value'', 0);' ...
                    'set(findobj(''parent'', gcbf, ''tag'', ''spec_no''), ''value'', 0);'  ...
                    'set(findobj(''parent'', gcbf, ''tag'', ''spec_noall''), ''value'', 0);'];
                 set_no =  [ 'set(findobj(''parent'', gcbf, ''tag'', ''spec_yesall''), ''value'', 0);' ...
                    'set(findobj(''parent'', gcbf, ''tag'', ''spec_yes''), ''value'', 0);' ...
                    'set(findobj(''parent'', gcbf, ''tag'', ''spec_noall''), ''value'', 0);'];
                 set_noall =  [ 'set(findobj(''parent'', gcbf, ''tag'', ''spec_yesall''), ''value'', 0);' ...
                    'set(findobj(''parent'', gcbf, ''tag'', ''spec_yes''), ''value'', 0);' ...
                    'set(findobj(''parent'', gcbf, ''tag'', ''spec_no''), ''value'', 0);'  ];
                 spec_ans = inputgui({[1] [1] [1 1 ] [1 1] [1]}, ...
                         { {'style' 'text' 'string' ['Spectrum infomation [' num2str(round(fave(1))) ' ' num2str(round(fave(end)))  '] Hz already exists for dataset: '  EEG.setname '. ' ] } ...
                         {'style' 'text' 'string' 'Would you like to recalculate the spectrum and overwrite these information?' } ...
                         {'style' 'checkbox' 'tag' 'spec_yes' 'string' 'Yes' 'value' 1 'Callback' set_yes }  ...
                         {'style' 'checkbox' 'tag' 'spec_yesall' 'string' 'Yes to all datasets' 'value' 0 'Callback' set_yesall }  ...
                         {'style' 'checkbox' 'tag' 'spec_no' 'string' 'Use existing spectrum info' 'value' 0 'Callback' set_no } ...
                         {'style' 'checkbox' 'tag' 'spec_noall' 'string' 'Use existing info for all sets' 'value' 0 'Callback' set_noall } {} }, ...
                     '', 'Recalculate spectrum information -- part of cls_spec()'); 
                 switch (find(celltomat(spec_ans)))
                     case 1
                         [EEG_etc, X, f, overwrite] = cls_spec(EEG, comp, freqrange, arg, 1); % overwrite the info in this dataset
                         overwrite = 0; %but ask before overwriting the rest of the datasets
                         return;
                     case 2
                         [EEG_etc, X, f, overwrite] = cls_spec(EEG, comp, freqrange, arg, 1);
                         return;
                     case 3
                         overwrite = 0; % don't overwrite the info in this datase but keep asking for the other datasets
                     case 4
                         overwrite = 2; % keep the info in all datasets
                 end
             end
             
             X = zeros(length(comp),md) ;
             for k = 1:length(comp)
                 X(k,:) = floatread([ EEG.etc.icaspec 'm'], [md 1],[],md*comp(k))';
             end
             f = fave;
             eval ( ['cd '  olddir]);
             return
         else % overwrite existing spectrum (use exisiting spectrum but take new frequency boundaries). 
             f = floatread(EEG.etc.icaspec, [d 1]);
             if ~isempty(freqrange)
                 maxind = max(find(f <= freqrange(end)));
                 minind = min(find(f >= freqrange(1)));
                 if f(end) < freqrange(end)
                     disp(['Warning! Requested high frequency limit, ' num2str(freqrange(end)) 'Hz, is out of bound.']);
                 end
                 if f(1) > freqrange(1)
                     disp(['Warning! Requested low frequency limit, ' num2str(freqrange(1)) 'Hz, is out of bound.']);
                 end
                 f = f(minind:maxind);
             else
                 minind = 1;
                 maxind = d;
             end
             X = zeros(length(comp),maxind-minind+1) ;
             for k = 1:length(comp)
                 tmp = floatread(EEG.etc.icaspec, [d 1],[],d*comp(k));
                 X(k,:) =  tmp(minind:maxind)';
             end
             X = X - mean(X,2)*ones(1,length(f)); %remove mean
             X = X - ones(size(X,1),1)*mean(X); %remove the mean from each frequency across all components
             if minind ~= 1 | maxind ~= d %new removed mean values
                 if ~isempty(arg)
                    EEG.etc.icaspecmparams = {length(f), arg{:}};
				else
                    EEG.etc.icaspecmparams = {length(f)};
				end
                %Save the updated dataset
				try
                    EEG = pop_saveset( EEG, 'filename', EEG.filename, 'filepath', EEG.filepath, 'savemode','twofiles');
				catch,
                    error([ 'cls_spec: problems saving into path ' EEG.filepath])
				end
                EEG_etc = EEG.etc;
                floatwrite([f X'], [ EEG.filename(1:end-3) 'icaspecm']);%save removed mean spectra info in float file
             end
             eval ([ 'cd ' olddir]); 
             return
         end
     end
 end
 
%no spectrum information
if isempty(EEG.icaact)
    EEG = eeg_checkset( EEG, 'loaddata' ); %load EEG.data and EEG.icaact
end
if ~isempty(arg)
    [X, f] = spectopo(EEG.icaact, EEG.pnts, EEG.srate, 'plot', 'off', arg);  
else
    [X, f] = spectopo(EEG.icaact, EEG.pnts, EEG.srate, 'plot', 'off');
end

%save spectrum in file
olddir = pwd;
eval ( ['cd '  EEG.filepath]);
floatwrite([f X'], [ EEG.filename(1:end-3) 'icaspec']); %save spectra info in float file

%update the info in the dataset
EEG.etc.icaspec = [ EEG.filename(1:end-3) 'icaspec'];
if ~isempty(arg)
    EEG.etc.icaspecparams = {length(f), arg{:}};
else
    EEG.etc.icaspecparams = {length(f)};
end

%Slect desired components
X = X(comp,:); 
if ~isempty('freqrange')
    maxind = max(find(f <= freqrange(end)));
    minind = min(find(f >= freqrange(1)));
    if f(end) < freqrange(end)
        disp(['Warning! Requested high frequency limit, ' num2str(freqrange(end)) 'Hz, is out of bound.']);
    end
    if f(1) > freqrange(1)
        disp(['Warning! Requested low frequency limit, ' num2str(freqrange(1)) 'Hz, is out of bound.']);
    end
    f = f(minind:maxind);
    X = X(:,minind:maxind);
end               
X = X - mean(X,2)*ones(1,length(f)); %remove mean
X = X - ones(size(X,1),1)*mean(X); %remove the mean from each frequency across all components
if ~isempty(arg)
    EEG.etc.icaspecmparams = {length(f), arg{:}};
else
    EEG.etc.icaspecmparams = {length(f)};
end
floatwrite([f X'], [ EEG.filename(1:end-3) 'icaspecm']);%save removed mean spectra info in float file

%Save the updated dataset
try
    EEG = pop_saveset( EEG, 'savemode', 'resave');
catch,
    error([ 'cls_spec: problems saving into path ' EEG.filepath])
end
EEG_etc = EEG.etc;
eval ([ 'cd ' olddir]); 
