% biosig2eeglab() - convert BIOSIG structure to EEGLAB structure
%
% Usage:
%   >> OUTEEG = pop_biosig2eeglab(hdr, data, interval);
%
% Inputs:
%   hdr   - BIOSIG header
%   data  - BIOSIG data array
%
% Optional input:
%   interval - BIOSIG does not remove event which are outside of
%              the data range when importing data range subsets. This 
%              parameter helps fix this problem.
%
% Outputs:
%   OUTEEG   - EEGLAB data structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, Oct. 29, 2009-
%
% Note: BIOSIG toolbox must be installed. Download BIOSIG at 
%       http://biosig.sourceforge.net
%       Contact a.schloegl@ieee.org for troubleshooting using BIOSIG.

% Copyright (C) 2003 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
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

function EEG = biosig2eeglab(dat, DAT, interval, channels, importevent, importannot)

if nargin < 2
    help biosig2eeglab;
    return;
end
if nargin < 3
    interval = [];
end
if nargin < 4
    channels = [];
end
if nargin < 5
    importevent = 0;
end
if nargin < 6
    importannot = 0;
end

% import data
% -----------
if abs(dat.NS - size(DAT,1)) > 1
    DAT = DAT';
end
EEG = eeg_emptyset;

% convert to seconds for sread
% ----------------------------
if max(dat.InChanSelect) > size(DAT,1)
    dat.InChanSelect = [1:size(DAT,1)];
end
EEG.nbchan          = length(dat.InChanSelect); %= size(DAT,1);
EEG.srate           = dat.SampleRate(1);
EEG.data            = DAT(dat.InChanSelect,:);  %DAT
clear DAT;
% try  % why would you do the following???????  JO
%     EEG.data            = EEG.data';
% catch,
%     pack;
%     EEG.data            = EEG.data';
% end;    
EEG.setname 		= sprintf('%s file', dat.TYPE);
EEG.comments        = [ 'Original file: ' dat.FileName ];
EEG.xmin            = 0;
nepoch              = dat.NRec;
EEG.trials   = nepoch;
EEG.pnts     = size(EEG.data,2)/nepoch;

if isfield(dat,'T0')
    EEG.etc.T0 = dat.T0; % added sjo
end

if isfield(dat, 'Label') && ~isempty(dat.Label)
    if ischar(dat.Label)
        EEG.chanlocs = struct('labels', cellstr(char(dat.Label(dat.InChanSelect)))); % 5/8/2104 insert (dat.InChanSelect) Ramon
    else
        % EEG.chanlocs = struct('labels', dat.Label(1:min(length(dat.Label), size(EEG.data,1))));
        EEG.chanlocs = struct('labels', dat.Label(dat.InChanSelect)); % sjo added 120907 to avoid error below % 5/8/2104 insert (dat.InChanSelect) Ramon
    end
    if length(EEG.chanlocs) > EEG.nbchan, EEG.chanlocs = EEG.chanlocs(1:EEG.nbchan); end
    
    % check BIOSEMI channel types
    if strcmpi(dat.TYPE, 'BDF')
        bdftype = { 'EXG1'  'EXG2'  'EXG3'  'EXG4'  'EXG5'  'EXG6'  'EXG7'  'EXG8'  'GSR1'  'GSR2'  'Erg1'  'Erg2'  'Resp'  'Plet'  'Temp';
                    'MISC'  'MISC'  'MISC'  'MISC'  'MISC'  'MISC'  'MISC'  'MISC'  'GSR'   'GSR'   'MISC'  'MISC'  'RESP'  'PPG'   'TEMP' };
        for iChan = 1:length(EEG.chanlocs)
            matchChan = strmatch(EEG.chanlocs(iChan).labels, bdftype(1,:), 'exact');
            if ~isempty(matchChan)
                EEG.chanlocs(iChan).type = bdftype{2,matchChan};
            else
                EEG.chanlocs(iChan).type = 'EEG';
            end
        end
    end
end
EEG = eeg_checkset(EEG);

% extract events % this part I totally revamped to work...  JO
% --------------
EEG.event = [];

% startval = mode(EEG.data(end,:)); % my code
% for p = 2:size(EEG.data,2)-1
%     [codeout] = code(EEG.data(end,p));
%     if EEG.data(end,p) > EEG.data(end,p-1) & EEG.data(end,p) >= EEG.data(end,p+1)
%         EEG.event(end+1).latency =  p;
%         EEG.event(end).type = bitand(double(EEG.data(end,p)-startval),255);
%     end
% end

% lastout = mod(EEG.data(end,1),256);newevs = []; % andrey's code 8 bits
% codeout = mod(EEG.data(end,2),256);
% for p = 2:size(EEG.data,2)-1
%     nextcode = mod(EEG.data(end,p+1),256);
%     if codeout > lastout & codeout >= nextcode
%         newevs = [newevs codeout];
%         EEG.event(end+1).latency =  p;
%         EEG.event(end).type = codeout;
%     end
%     lastout = codeout;
%     codeout = nextcode;
% end

%lastout = mod(EEG.data(end,1),256*256);newevs = []; % andrey's code 16 bits
%codeout = mod(EEG.data(end,2),256*256);
%for p = 2:size(EEG.data,2)-1
%    nextcode = mod(EEG.data(end,p+1),256*256);
%    if (codeout > lastout) & (codeout >= nextcode)
%        newevs = [newevs codeout];
%        EEG.event(end+1).latency =  p;
%        EEG.event(end).type = codeout;
%    end
%    lastout = codeout;
%    codeout = nextcode;
%end

% if strcmp(dat.TYPE,'EDF') % sjo added 120907
%     disp('filetype EDF does not support events');
%     importevent = 0;
% end

if importevent
    if isfield(dat, 'BDF')
        if dat.BDF.Status.Channel <= size(EEG.data,1)
            EEG.data(dat.BDF.Status.Channel,:) = [];
        end
        EEG.nbchan = size(EEG.data,1);
        if ~isempty(EEG.chanlocs) && dat.BDF.Status.Channel <= length(EEG.chanlocs)
            EEG.chanlocs(dat.BDF.Status.Channel,:) = [];
        end
    elseif isempty(dat.EVENT.POS)
        disp('Extracting events from last EEG channel...');
        %Modifieded by Andrey (Aug.5,2008) to detect all non-zero codes: 
        if length(unique(EEG.data(end, 1:100))) > 20
            disp('Warning: event extraction failure, the last channel contains data');
        elseif length(unique(EEG.data(end, :))) > 1000
            disp('Warning: event extraction failure, the last channel contains data');
        else
            thiscode = 0;
            for p = 1:size(EEG.data,2)*size(EEG.data,3)-1
                prevcode = thiscode;
                thiscode = mod(EEG.data(end,p),256*256);   % andrey's code - 16 bits 
                if (thiscode ~= 0) && (thiscode~=prevcode) 
                    EEG.event(end+1).latency =  p;
                    EEG.event(end).type = thiscode;
                end
            end
            EEG.data(end,:) = [];
            EEG.chanlocs(end) = [];
        end
        % recreate the epoch field if necessary
        % -------------------------------------
        if EEG.trials > 1
            for i = 1:length(EEG.event)
                EEG.event(i).epoch = ceil(EEG.event(i).latency/EEG.pnts);
            end
        end
        EEG = eeg_checkset(EEG, 'eventconsistency');
    end

    if ~isempty(dat.EVENT.POS)    
        if isfield(dat, 'out') % Alois fix for event interval does not work
            if isfield(dat.out, 'EVENT')
                dat.EVENT = dat.out.EVENT;
            end
        end
        EEG.event = biosig2eeglabevent(dat.EVENT, interval, importannot); % Toby's fix

        % recreate the epoch field if necessary
        % -------------------------------------
        if EEG.trials > 1
            for i = 1:length(EEG.event)
                EEG.event(i).epoch = ceil(EEG.event(i).latency/EEG.pnts);
            end
        end

        EEG = eeg_checkset(EEG, 'eventconsistency');
    elseif isempty(EEG.event) 
        disp('Warning: no event found. Events might be embedded in a data channel.');
        disp('         To extract events, use menu File > Import Event Info > From data channel');
    end
end

% convert data to single if necessary
% -----------------------------------
EEG = eeg_checkset(EEG,'makeur');   % Make EEG.urevent field
EEG = eeg_checkset(EEG);
