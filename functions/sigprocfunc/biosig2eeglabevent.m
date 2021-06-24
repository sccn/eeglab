% biosig2eeglabevent() - convert biosig events to EEGLAB event structure
%
% Usage:
%   >> eeglabevent = biosig2eeglabevent( biosigevent, interval)
%
% Inputs:
%   biosigevent    - BioSig event structure
%   interval       - Period to extract events for, in frames.
%                    Default [] is all.
%
% Outputs:
%   eeglabevent    - EEGLAB event structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2006-

% Copyright (C) 13 2006- Arnaud Delorme, Salk Institute, arno@salk.edu
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

function event = biosig2eeglabevent(EVENT, interval, importEDFplus)

if nargin < 2
    interval = [];
end
if nargin < 3
    importEDFplus = false;
end

event = [];
disp('Importing data events...');

EVT=sopen('eventcodes.txt');sclose(EVT);

% If the interval variable is empty, import all events.
if isempty(interval)
    if isfield(EVENT, 'Teeg')
        event = EVENT.Teeg;
    end
    if isfield(EVENT, 'TYP')
        for index = 1:length( EVENT.TYP )
            eType = EVENT.TYP(index);
            
            % use file in
            % https://sccn.ucsd.edu/bugzilla/show_bug.cgi?id=1387 to test
            % for boundary events
            if eType < 256 && importEDFplus && isfield(EVENT,'CodeDesc') && eType < length(EVENT.CodeDesc)
                event(index).type = EVENT.CodeDesc{eType};
                event(index).edftype = eType;
            elseif isfield(EVT, 'EVENT') && isfield(EVT.EVENT,'CodeIndex') && isfield(EVT.EVENT,'CodeDesc') && importEDFplus
                try
                    event(index).type = EVT.EVENT.CodeDesc{EVT.EVENT.CodeIndex==eType};
                    event(index).edftype = eType;
                catch
                    event(index).type = eType;
                end
                if eType == 32766 || eType == 32767
                    event(index).edfplustype = event(index).type;
                    event(index).type = 'boundary';
                end
            else
                event(index).type = eType;
            end
        end
    end
    if isfield(EVENT, 'POS')
        for index = 1:length( EVENT.POS )
            event(index).latency = EVENT.POS(index);
        end
    end
    if isfield(EVENT, 'DUR')
        if any( [ EVENT.DUR ] )
            for index = 1:length( EVENT.DUR )
                event(index).duration = EVENT.DUR(index);
            end
        end
    end
    if isfield(EVENT, 'CHN')
        if any( [ EVENT.CHN ] )
            for index = 1:length( EVENT.CHN )
                event(index).chanindex = EVENT.CHN(index);
            end
        end
    end
% If a subinterval was specified, select only events that fall in that range, and
% edit duration field if it exceeds that range.
elseif isfield(EVENT,'POS')
    count = 1;
    for index = 1:length(EVENT.POS)
        pos_tmp = EVENT.POS(index) - interval(1) + 1;
        if pos_tmp > 0 && EVENT.POS(index) <= interval(2)
            event(count).latency = pos_tmp;
            if isfield(EVENT, 'TYP')
                eType = EVENT.TYP(index);

                if isfield(EVENT, 'CodeDesc')
                    if eType < 256 && importEDFplus && eType < length(EVENT.CodeDesc)
                        event(index).type = EVENT.CodeDesc{eType};
                        event(index).edftype = eType;
                    elseif isfield(EVT, 'EVENT') && isfield(EVT.EVENT,'CodeIndex') && isfield(EVT.EVENT,'CodeDesc') && importEDFplus
                        event(index).type = EVT.EVENT.CodeDesc{EVT.EVENT.CodeIndex==eType};
                        event(index).edftype = eType;
                        if eType == 32766 || eType == 32767
                            event(index).edfplustype = event(index).type;
                            event(index).type = 'boundary';
                        end
                    else
                        event(index).type = eType;
                    end
                else
                    event(index).type = eType;
                end
            end
            if isfield(EVENT, 'CHN')
                event(count).chanindex = EVENT.CHN(index);
            end
            if isfield(EVENT, 'DUR')
                event(count).duration = min(EVENT.DUR(index), interval(2) - EVENT.POS(index));
            end
            count = count + 1;
        end
    end
end
