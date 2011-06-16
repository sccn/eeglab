% eeg_checkchanlocs() - Check the consistency of the channel locations structure 
%                  of an EEGLAB dataset.
%
% Usage:
%  >> EEG = eeg_checkchanlocs( EEG, 'key1', value1, 'key2', value2, ... ); 
%  >> [chanlocs chaninfo] = eeg_checkchanlocs( chanlocs, chaninfo, 'key1', value1, 'key2', value2, ... ); 
%
% Inputs:
%   EEG      - EEG dataset
%   chanlocs - EEG.chanlocs structure
%   chaninfo - EEG.chaninfo structure
%
% Outputs:
%   EEG        - new EEGLAB dataset with updated channel location structures 
%                EEG.chanlocs, EEG.urchanlocs, EEG.chaninfo
%   chanlocs   - updated channel location structure
%   chaninfo   - updated chaninfo structure
%
% Author: Arnaud Delorme, SCCN/INC/UCSD, March 2, 2011

% Copyright (C) SCCN/INC/UCSD, March 2, 2011, arno@salk.edu
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

function [chans, chaninfo, chanedit]= eeg_checkchanlocs(chans, chaninfo);

if nargin < 1 
    help eeg_checkchanlocs;
    return;
end;

if nargin < 2
    chaninfo = [];
end;

processingEEGstruct = 0;
if isfield(chans, 'data')
    processingEEGstruct = 1;
    tmpEEG = chans;
    chans = tmpEEG.chanlocs;
    chaninfo = tmpEEG.chaninfo;
end;

if ~isfield(chans, 'datachan')
    chanedit = insertchans(chans, chaninfo);
else
    chanedit = chans;
end;
nosevals       = { '+X' '-X' '+Y' '-Y' };
if ~isfield(chaninfo, 'plotrad'), chaninfo.plotrad = []; end;
if ~isfield(chaninfo, 'shrink'),  chaninfo.shrink = [];  end;
if ~isfield(chaninfo, 'nosedir'), chaninfo.nosedir = nosevals{1}; end;

% handles deprecated fields
% -------------------------
plotrad  = [];
if isfield(chanedit, 'plotrad'),
    plotrad = chanedit(1).plotrad;
    chanedit = rmfield(chanedit, 'plotrad');
    if isstr(plotrad) & ~isempty(str2num(plotrad)), plotrad = str2num(plotrad); end;
    chaninfo.plotrad = plotrad;
end;
if isfield(chanedit, 'shrink') && ~isempty(chanedit(1).shrink)
    shrinkorskirt = 1;
    if ~isstr(chanedit(1).shrink)
        plotrad = 0.5/(1-chanedit(1).shrink); % convert old values
    end;
    chanedit = rmfield(chanedit, 'shrink');
    chaninfo.plotrad = plotrad;
end;

% set non-existant fields to []
% -----------------------------
fields    = { 'labels' 'theta' 'radius' 'X'   'Y'   'Z'   'sph_theta' 'sph_phi' 'sph_radius' 'type' 'ref' 'urchan' };
fieldtype = { 'str'    'num'   'num'    'num' 'num' 'num' 'num'       'num'     'num'        'str'  'str' 'num'    };
if ~isempty(chanedit)
    for index = 1:length(fields)
        if ~isfield(chanedit, fields{index})
            % new field
            % ---------
            if strcmpi(fieldtype{index}, 'num')
                chanedit = setfield(chanedit, {1}, fields{index}, []);
            else
                for indchan = 1:length(chanedit)
                    chanedit = setfield(chanedit, {indchan}, fields{index}, '');
                end;
            end;
        else
            % existing fields
            % ---------------
            eval([ 'allvals = { chanedit.' fields{index} '};' ] );
            if strcmpi(fieldtype{index}, 'num')
                numok = cellfun(@isnumeric, allvals);
                if any(numok == 0)
                    for indConvert = find(numok == 0)
                        chanedit = setfield(chanedit, {indConvert}, fields{index}, []);
                    end;
                end;
            else
                strok = cellfun(@isstr, allvals);
                if strcmpi(fields{index}, 'labels'), prefix = 'E'; else prefix = ''; end;
                if any(strok == 0)
                    for indConvert = find(strok == 0)
                        try
                              strval   = [ prefix num2str(getfield(chanedit, {indConvert}, fields{index})) ];
                              chanedit = setfield(chanedit, {indConvert}, fields{index}, strval);
                        catch chanedit = setfield(chanedit, {indConvert}, fields{index}, '');
                        end;
                    end;
                end;
            end;
        end;
    end;
end;
if exist('orderfields') == 2
    try,
        chanedit = orderfields(chanedit, fields);
    catch, end;
end;

% check if duplicate channel label
% --------------------------------
if isfield(chanedit, 'labels')
    if length( { chanedit.labels } ) > length( unique({ chanedit.labels } ) )
        disp('Warning: some channels have the same label');
    end;
end;

% remove fields
% -------------
if isfield(chanedit, 'sph_phi_besa'  ), chanedit = rmfield(chanedit, 'sph_phi_besa'); end;
if isfield(chanedit, 'sph_theta_besa'), chanedit = rmfield(chanedit, 'sph_theta_besa'); end;

% reconstruct the chans structure
% -------------------------------
[chans chaninfo.nodatchans] = getnodatchan( chanedit );
if ~isfield(chaninfo, 'nodatchans'), chaninfo.nodatchans = []; end;
if isempty(chanedit)
    for iField = 1:length(fields)
        chanedit = setfield(chanedit, fields{iField}, []);
    end;
end;

if processingEEGstruct
    tmpEEG.chanlocs = chans;
    tmpEEG.chaninfo = chaninfo;
    chans = tmpEEG;
end;

% ---------------------------------------------
% separate data channels from non-data channels
% ---------------------------------------------
function [chans, fids] = getnodatchan(chans)
if isfield(chans, 'datachan')
    for ind = 1:length(chans)
        if isempty(chans(ind).datachan)
            chans(ind).datachan = 0;
        end;
    end;
    alldatchans = [ chans.datachan ];
    chans = rmfield(chans, 'datachan');
    fids  = chans(find(alldatchans == 0));
    chans = chans(find(alldatchans));
else
    fids = [];
end;

% ----------------------------------------
% fuse data channels and non-data channels
% ----------------------------------------
function [chans, chaninfo] = insertchans(chans, chaninfo, nchans)
if nargin < 3, nchans = length(chans); end;
[chans.datachan] = deal(1);
if isfield(chans,'type')
    mask = strcmpi({chans.type},'FID') | strcmpi({chans.type},'IGNORE');
    [chans(mask).datachan] = deal(0);
end
if length(chans) > nchans & nchans ~= 0 % reference at the end of the structure
    chans(end).datachan = 0;
end;
if isfield(chaninfo, 'nodatchans')
    if ~isempty(chaninfo.nodatchans)
        chanlen = length(chans);
        for index = 1:length(chaninfo.nodatchans)
            fields = fieldnames( chaninfo.nodatchans );
            ind = chanlen+index;
            for f = 1:length( fields )
                chans = setfield(chans, { ind }, fields{f}, getfield( chaninfo.nodatchans, { index },  fields{f}));
            end;
            chans(ind).datachan = 0;
        end;
        chaninfo = rmfield(chaninfo, 'nodatchans');
        
        % put these channels first
        % ------------------------
        % tmp = chans(chanlen+1:end);
        % chans(length(tmp)+1:end) = chans(1:end-length(tmp));
        % chans(1:length(tmp)) = tmp;
    end;
end;
