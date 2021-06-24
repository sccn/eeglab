% eeg_getversion() - obtain EEGLAB version number (version is embeded in
%                    the script, edit the function to see the version).
%
% Usage:
%     >> vers = eeg_getversion;
%     >> [vers vnum] = eeg_getversion;
%
% Outputs:
%    vers = [string] EEGLAB version number
%    vnum = [float] numerical value for the version. For example 11.3.2.4b
%           is converted to 11.324
%
% Authors: Arnaud Delorme, SCCN/INC/UCSD, 2010

% Copyright (C) 2010  Arnaud Delorme, SCCN/INC/UCSD
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

function [vers, versnum, releaseDate] = eeg_getversion

vers        = '2021.0';
releaseDate = ''; % 30-Apr-19 14:55:42; unix date -> date +"%d-%b-%y %T"

% get numerical version number
tmpvers = vers;
if isnan(str2double(tmpvers(end))), tmpvers(end) = []; end
indsDot = find(tmpvers == '.' );
tmpvers(indsDot(2:end)) = [];
versnum = str2double(tmpvers);
if isnan(versnum), versnum = []; end
