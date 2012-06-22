% std_erpimage() - Compute ERP images and save them on disk.
%
% Usage:
%   >> std_erpimage( EEG, 'key', 'val', ...);
%
% Inputs:
%   EEG          - a loaded epoched EEG dataset structure. May be an array
%                  of such structure containing several datasets.
%
% Optional inputs:
%   'components' - [numeric vector] components of the EEG structure for which 
%                  the measure will be computed {default|[] -> all}
%   'channels'   - [cell array] channels of the EEG structure for which 
%                  activation ERPs will be computed {default|[] -> none}
%   'trialindices' - [cell array] indices of trials for each dataset.
%                  Default is all trials.
%   'recompute'  - ['on'|'off'] force recomputing data file even if it is 
%                  already on disk.
%   'rmcomps'    - [integer array] remove artifactual components (this entry
%                  is ignored when plotting components). This entry contains 
%                  the indices of the components to be removed. Default is none.
%   'interp'     - [struct] channel location structure containing electrode
%                  to interpolate ((this entry is ignored when plotting 
%                  components). Default is no interpolation.
%   'fileout'    - [string] name of the file to save on disk. The default
%                  is the same name (with a different extension) as the 
%                  dataset given as input.
%
% ERPimage options:
%   'concatenate' - ['on'|'off'] concatenate single trial of different
%                  subjects for plotting ERPimages ('on'). The default
%                  ('off') computes an ERPimage for each subject and then
%                  averages these ERPimages. This allows to perform
%                  statistics (the 'on' options does not allow statistics).
%   'smoothing'  - Smoothing parameter (number of trials). {Default: 10}
%                  erpimage() equivalent: 'avewidth'
%   'nlines'     - Number of lines for ERPimage. erpaimge() equivalent is 
%                  'decimate'. Note that this parameter must be larger than
%                  the minimum number of trials in each design cell 
%                  {Default: 10}
%   'sorttype'   - Sorting event type(s) ([int vector]; []=all). See Notes below.
%                  Either a string or an integer.
%   'sortwin'    - Sorting event window [start, end] in seconds ([]=whole epoch)
%   'sortfield'  - Sorting field name. {default: none}. See Notes below.
%   'erpimageopt'  - erpimage() options, separated by commas (Ex: 'erp', 'cbar').
%                  {Default: none}. For further details see >> erpimage help
% Outputs:
%   erpimagestruct - structure containing ERPimage information that is
%                    been saved on disk.
%
%   Files are saved on disk.
%    [dataset_file].icaerpim     % component ERPimage file
% OR
%    [dataset_file].daterpim     % channel ERPimage file
%
% Author: Arnaud Delorme, SCCN & CERCO, CNRS, 2011-

% Copyright (C) 2011 Arnaud Delorme
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

function allerpimage = std_erpimage( EEG, varargin);

if nargin < 1
    help std_erpimage;
    return;
end

allerpimage = [];
[opt moreopts] = finputcheck( varargin, { ...
    'components'    'integer'     []                    [];
    'channels'      { 'cell','integer' }  { [] [] }     {};
    'trialindices' { 'integer','cell' }   []            [];
    'recompute'     'string'      { 'on','off' }        'off';
    'savefile'      'string'      { 'on','off' }        'on';
    'fileout'       'string'      []                    '';
    'rmcomps'       'cell'        []                    cell(1,length(EEG));
    'interp'        'struct'      { }                   struct([]);
    'nlines'        ''            []                    10;
    'smoothing'     ''            []                    10;
    'sorttype'       ''           {}                    '';
    'sortwin'        ''           {}                    [];
    'sortfield'      ''           {}                    '';
    'concatenate'   'string'      { 'on','off' }        'off';
    'erpimageopt'   'cell'        {}                    {}}, ...
    'std_erpimage', 'ignore');
if isstr(opt), error(opt); end;
if length(EEG) == 1 && isempty(opt.trialindices), opt.trialindices = { [1:EEG.trials] }; end;
if isempty(opt.trialindices), opt.trialindices = cell(length(EEG)); end;
if ~iscell(opt.trialindices), opt.trialindices = { opt.trialindices }; end;
if isfield(EEG,'icaweights')
    numc = size(EEG(1).icaweights,1);
else
    error('EEG.icaweights not found');
end
if isempty(opt.components)
    opt.components = 1:numc;
end

% filename
% --------
if isempty(opt.fileout), opt.fileout = fullfile(EEG(1).filepath, EEG(1).filename(1:end-4)); end;
if ~isempty(opt.channels)
    filenameshort = [ opt.fileout '.daterpim'];
    prefix = 'chan';
    if iscell(opt.channels)
        if ~isempty(opt.interp)
            opt.indices = eeg_chaninds(opt.interp, opt.channels, 0);
        else
            opt.indices = eeg_chaninds(EEG(1), opt.channels, 0);
            for ind = 2:length(EEG)
                if ~isequal(eeg_chaninds(EEG(ind), opt.channels, 0), opt.indices)
                    error([ 'Channel information must be consistant when ' 10 'several datasets are merged for a specific design' ]);
                end;
            end;
        end;
    else
        opt.indices = opt.channels;
    end;
else
    opt.indices = opt.components;
    filenameshort = [ opt.fileout '.icaerpim'];
    prefix = 'comp';
end;
filename = filenameshort;

% ERP information found in datasets
% ---------------------------------
if exist(filename) && strcmpi(opt.recompute, 'off')
    fprintf('File "%s" found on disk, no need to recompute\n', filenameshort);
    return;
end

allerpimage = [];
if strcmpi(opt.concatenate, 'off')
    % compute ERP images
    % ------------------
    if isempty(opt.channels)
         X = eeg_getdatact(EEG, 'component', opt.indices, 'trialindices', opt.trialindices );
    else X = eeg_getdatact(EEG, 'channel'  , opt.indices, 'trialindices', opt.trialindices, 'rmcomps', opt.rmcomps, 'interp', opt.interp);
    end;
    if ~isempty(opt.sorttype)
         events = eeg_getepochevent(EEG, 'type', opt.sorttype, 'timewin', opt.sortwin, 'fieldname', opt.sortfield, 'trials', opt.trialindices);
    else events = [];
    end;
        
    % reverse engeeneering the number of lines for ERPimage
    finallines = opt.nlines;
    if ~isempty(events)
         lastx  = sum(~isnan(events));
    else lastx  = size(X,3);
    end;
    if lastx < finallines + floor((opt.smoothing-1)/2) + 3
        error('The default number of ERPimage lines is too large for one of the dataset');
    end;
    firstx = 1;
    xwidth = opt.smoothing;
    %xadv   = lastx/finallines;
    nout   = finallines; %floor(((lastx-firstx+xadv+1)-xwidth)/xadv);
    nlines = (lastx-xwidth)/(nout-0.5)*i; % make it imaginary
    %nlines = ceil(lastx/((lastx-firstx+1-xwidth)/(nout-1)));


    if 0
        % testing conversion back and forth
        % ---------------------------------
        for lastx = 20:300
            for xwidth = 1:19
                for nlines = (xwidth+1):100
                    
                    nout = floor(((lastx+nlines)-xwidth)/nlines);
                    realnlines = (lastx-xwidth)/(nout-0.5);
                    noutreal = floor(((lastx+realnlines)-xwidth)/realnlines);
                    
                    if nout ~= noutreal
                        error('Wrong conversion 2');
                    end;
                    
                end;
            end;
        end;
    end;
    
    for index = 1:size(X,1)
        [tmpX tmpevents] = erpimage(squeeze(X(index,:,:)), events, EEG(1).times, '', opt.smoothing, nlines, 'noplot', 'on', opt.erpimageopt{:}, moreopts{:});
        allerpimage = setfield(allerpimage, 'events', tmpevents);
        allerpimage = setfield(allerpimage, [ prefix int2str(opt.indices(index)) ], tmpX');
    end;
else
    % generate dynamic loading commands
    % ---------------------------------
    for dat = 1:length(EEG)
        filenames{dat} = fullfile(EEG(1).filepath, EEG(1).filename);
    end;
    allerpimage.times = EEG(1).times;
    for index = 1:length(opt.indices)
        if ~isempty(opt.channels)
             com = sprintf('squeeze(eeg_getdatact(%s, ''interp'', chanlocsforinterp));', vararg2str( { filenames 'channel'  , opt.indices(index), 'rmcomps', opt.rmcomps, 'trialindices', opt.trialindices } ));
        else com = sprintf('squeeze(eeg_getdatact(%s));', vararg2str( { filenames 'component', opt.indices(index), 'trialindices', opt.trialindices } ));
        end;
        allerpimage = setfield(allerpimage, [ prefix int2str(opt.indices(index)) ], com);
    end;
    allerpimage = setfield(allerpimage, 'chanlocsforinterp', opt.interp);
    if ~isempty(opt.sorttype)
         events = eeg_getepochevent(EEG, 'type', opt.sorttype, 'timewin', opt.sortwin, 'fieldname', opt.sortfield, 'trials', opt.trialindices);
         %geteventcom = sprintf('eeg_getepochevent(%s);', vararg2str( { filenames 'type', opt.sorttype, 'timewin', opt.sortwin, 'fieldname', opt.sortfield } ));
    else events = [];
    end;
    allerpimage = setfield(allerpimage, 'events', events);
end;
allerpimage.times       = EEG(1).times;
allerpimage.parameters  = varargin;
allerpimage.datatype    = 'ERPIMAGE';
allerpimage.datafiles   = computeFullFileName( { EEG.filepath }, { EEG.filename });
allerpimage.datatrials  = opt.trialindices;

% Save ERPimages in file (all components or channels)
% ----------------------------------------------
if strcmpi(opt.savefile, 'on')
    if strcmpi(prefix, 'comp')
        std_savedat(filename, allerpimage);
    else
        tmpchanlocs = EEG(1).chanlocs;
        allerpimage.labels = opt.channels;
        std_savedat(filename, allerpimage);
    end;
end;

% compute full file names
% -----------------------
function res = computeFullFileName(filePaths, fileNames);
for index = 1:length(fileNames)
    res{index} = fullfile(filePaths{index}, fileNames{index});
end;
