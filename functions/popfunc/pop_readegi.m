% pop_readegi() - load a EGI EEG file (pop out window if no arguments).
%
% Usage:
%   >> EEG = pop_readegi;             % a window pops up
%   >> EEG = pop_readegi( filename );
%   >> EEG = pop_readegi( filename, datachunks, forceversion, fileloc);
%
% Inputs:
%   filename       - EGI file name
%   datachunks     - desired frame numbers (see readegi() help)
%                    option available from the command line only
%   forceversion   - [integer] force reading a specfic file version
%   fileloc        - [string] channel location file name. Default is
%                    'auto' (autodetection)
%
% Outputs:
%   EEG            - EEGLAB data structure
%
% Author: Arnaud Delorme, CNL / Salk Institute, 12 Nov 2002
%
% See also: eeglab(), readegi(), readegihdr()

% Copyright (C) 12 Nov 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [EEG, command] = pop_readegi(filename, datachunks, forceversion, fileloc);

EEG = [];
command = '';
if nargin < 4, fileloc = 'auto'; end
disp('Warning: This function can only import continuous files or');
disp('         epoch files with only one length for data epochs');

if nargin < 1
    % ask user
    [filename, filepath] = uigetfile('*.RAW;*.raw', ...
        'Choose an EGI RAW file -- pop_readegi()');
    drawnow;
    if filename == 0 return; end
    filename = [filepath filename];
    
    fid = fopen(filename, 'rb', 'b');
    if fid == -1, error('Cannot open file'); end
    head = readegihdr(fid); % read EGI file header
    fclose(fid);
    if head.segments ~= 0
        fileloc = '';
        floc = { 'GSN-HydroCel-32.sfp' 'GSN65v2_0.sfp' 'GSN129.sfp' 'GSN-HydroCel-257.sfp'};
        switch head.nchan
            case { 32 33 }, fileloc = floc;
            case { 64 65 }, fileloc =  {floc{2} floc{3:4} floc{1}};
            case { 128 129 }, fileloc = {floc{3} floc{4} floc{1:2}};
            case { 256 257 }, fileloc = {floc{4} floc{1:3}};
        end
        uilist = { { 'style' 'text' 'string' sprintf('Segment/frame number (default: 1:%d)', head.segments) } ...
                   { 'style' 'edit' 'string' '' } ...
                   { 'style' 'text' 'string' 'Channel location file (in eeglab/sample_locs)' } ...
                   { 'style' 'popupmenu' 'string' fileloc } ... 
                   { } ...
                   { 'style' 'text' 'string' ...
                   [ 'Note: Choosing the correct electrode location file for your data is critical.' 10 ...
        'Note that in some cases none of the channel location files listed will correspond' 10 ...
        'to your montage as EGI has different versions of caps and is creating new ones' 10 ...
        'constantly. Remember to check your montage in the channel editor and import. Channel' 10 ...
        'location files are stored in the sample_locs sub-folder of the EEGLAB distribution.' ] } };
        uigeometry = { [2 1] [2 1] [1] [1] };
        uigeomvert = [1 1 1.5 3];
        result = inputgui('uilist', uilist, 'geometry', uigeometry, 'geomvert', uigeomvert);
       
%         promptstr    = { sprintf('Segment/frame number (default: 1:%d)', head.segments) 'Channel location file (in eeglab/sample_locs)' };
%         inistr       = { '' fileloc(res{2})};
%         result       = inputdlg2( promptstr, 'Import EGI file -- pop_readegi()', 1,  inistr, 'pop_readegi');
        if length(result) == 0 return; end
        datachunks   = eval( [ '['  result{1} ']' ] );
        fileloc      = char(fileloc(result{2}));
    else
        datachunks   = [];
        disp('Only one segment, cannot read portion of the file');
    end
end

% load data
% ----------
EEG = eeg_emptyset;
if exist('datachunks') && exist('forceversion') && ~isempty(forceversion)
    [Head EEG.data Eventdata SegCatIndex] = readegi( filename, datachunks,forceversion);
elseif exist('forceversion') && ~isempty(forceversion)
    [Head EEG.data Eventdata SegCatIndex] = readegi( filename,[],forceversion);
elseif exist('datachunks')
    [Head EEG.data Eventdata SegCatIndex] = readegi( filename, datachunks);
    forceversion = [];
else
    [Head EEG.data Eventdata SegCatIndex] = readegi( filename);
    forceversion = [];
end
if ~isempty(Eventdata) && size(Eventdata,2) == size(EEG.data,2)
    EEG.data(end+1:end+size(Eventdata,1),:) = Eventdata;
end
EEG.comments        = [ 'Original file: ' filename ];
EEG.setname 		= 'EGI file';
EEG.nbchan          = size(EEG.data,1);
EEG.srate           = Head.samp_rate;
EEG.trials          = Head.segments;
EEG.pnts            = Head.segsamps;
EEG.xmin            = 0;

% importing the events
% --------------------
EEG = eeg_checkset(EEG);
if ~isempty(Eventdata)
    orinbchans = EEG.nbchan;
    for index = size(Eventdata,1):-1:1
        EEG = pop_chanevent( EEG, orinbchans-size(Eventdata,1)+index, 'edge', 'leading', ...
            'delevent', 'off', 'typename', Head.eventcode(index,:), ...
            'nbtype', 1, 'delchan', 'on');
        Head.eventcode(end,:) = [];
    end
    
    % renaming event codes
    % --------------------
    try,
        tmpevent = EEG.event;
        alltypes = { tmpevent.type };
        if ischar(alltypes{1})
            indepoc = strmatch('epoc', lower(alltypes), 'exact');
            indtim  = strmatch('tim0', lower(alltypes), 'exact');
            
            % if epoc but no tim0 then epoc represent pauses in recording
            if isempty(indtim) && ~isempty(indepoc)
                for index = indepoc
                    EEG.event(index).type = 'boundary';
                end
            end
            % other wise if both non-empty data epochs
            if ~isempty(indtim) && ~isempty(indepoc)
                if rem(size(EEG.data,2) / (length(indepoc)+1),1) == 0
                    EEG.event(index) = []; % remove epoch events
                    EEG.trials       = length(indepoc)+1;
                else
                    disp('Warning: data epochs detected but wrong data size');
                end
            end
        end
    catch, disp('Warning: event renaming failed'); end
end

% adding segment category indices
% -------------------------------
if ~isempty(SegCatIndex) && EEG.trials > 1
    try
        if ~isempty(EEG.event)
            for index = 1:length(EEG.event)
                EEG.event(index).category = Head.catname{SegCatIndex(EEG.event(index).epoch)};
            end
        else % create time-locking events
            for trial = 1:EEG.trials
                EEG.event(trial).epoch    = trial;
                EEG.event(trial).type     = 'TLE';
                EEG.event(trial).latency  = -EEG.xmin*EEG.srate+1+(trial-1)*EEG.pnts;
                EEG.event(trial).category = Head.catname{SegCatIndex(trial)};
            end
        end
    catch,
        disp('Warning: error while importing trial categories');
        EEG.event = rmfield(EEG.event, 'category');
    end
end
EEG = eeg_checkset(EEG, 'makeur');
EEG = eeg_checkset(EEG, 'eventconsistency');

% importing channel locations
% ---------------------------
if nargin < 1
    warndlg2( [ 'EEGLAB will now import a default electrode location file' 10 ...
        'for your data. Note that this might not correspond' 10 ...
        'to your montage as EGI has different versions of caps.' 10 ...
        'Check your montage in the channel editor and import' 10 ...
        'the correct location file if necessary.' ]);
end
if all(EEG.data(end,1:10) == 0)
    disp('Deleting empty data reference channel (reference channel location is retained)');
    EEG.data(end,:)   = [];
    EEG.nbchan        = size(EEG.data,1);
    EEG = eeg_checkset(EEG);
end
if ~isempty(fileloc)
    if strcmpi(fileloc, 'auto')
        EEG = readegilocs(EEG);
    else
        EEG = readegilocs(EEG, fileloc);
    end
end

if nargin < 1
    command = sprintf('EEG = pop_readegi(''%s'', %s);', filename, vararg2str({datachunks forceversion fileloc }) );
end

return;
