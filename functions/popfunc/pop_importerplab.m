% pop_importerplab() - import ERPLAB event list file and bin file into
%                      EEGLAB event structure for use in STUDY processing
%
% Usage:
%   >>  OUTEEG = pop_sample( INEEG, file1, file2);
%
% Inputs:
%   INEEG   - input EEG dataset
%   file1   - ERPLAB event list text file
%   file2   - ERPLAB bin file
%   ncbins  - [0|1] import all bins including non-contrast bins. Default is
%             0 (no).
%    
% Outputs:
%   OUTEEG  - output dataset

% Copyright (C) 2016 Arnaud Delorme
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

function [EEG, com] = pop_importerplab( EEG, file1, file2, ncbins)

% the command output is a hidden output that does not have to
% be described in the header

com = ''; % this initialization ensure that the function will return something
          % if the user press the cancel button            

% display help if not enough arguments
% ------------------------------------
if nargin < 1
	help pop_importerplab;
	return;
end
if nargin < 3
    ncbins = false;
end

% pop up window
% -------------
if nargin < 3
    promptstr    = { strvcat('Time-locking event type(s) ([]=all):', ...
                    'Select ''Edit > Event values'' to see type values.'), ...
                    'Epoch limits [start, end] in seconds:', ... 
                    'Name for the new dataset:', ... 
					'Out-of-bounds EEG rejection limits ([min max], []=none):'  };

    commandload1 = [ '[filename, filepath] = uigetfile(''*'', ''Select an event list text file'');' ...
        'if filename ~=0,' ...
        '   set(findobj(''parent'', gcbf, ''tag'', ''file1''), ''string'', [ filepath filename ]);' ...
        'end;' ...
        'clear filename filepath;' ];
    commandload2 = [ '[filename, filepath] = uigetfile(''*'', ''Select an bin text file'');' ...
        'if filename ~=0,' ...
        '   set(findobj(''parent'', gcbf, ''tag'', ''file2''), ''string'', [ filepath filename ]);' ...
        'end;' ...
        'clear filename filepath;' ];
   
    geometry = { [2 2 0.5] [2 2 0.5] [1] [1] };
    uilist = { { 'style' 'text'       'string' 'Select an event list text file ' } ...
              { 'style' 'edit'       'string' '' 'tag' 'file1' } ...
              { 'style' 'pushbutton' 'string' '...' 'callback' commandload1 } ...
              { 'style' 'text'       'string' 'Select a bin text file ' } ...
              { 'style' 'edit'       'string' '' 'tag' 'file2' } ...
              { 'style' 'pushbutton' 'string' '...' 'callback' commandload2 } ...
              {} ...
              { 'style' 'checkbox'   'string' 'Import all bins including non-contrast bins' 'tag' 'ncbins' } };
    [tmp1, tmp2, tmp3, results] = inputgui( geometry, uilist, 'pophelp(''pop_importerplab'');', 'Import ERPLAB info files -- pop_importerplab()' );

    file1  = results.file1;
    file2  = results.file2;
    ncbins = results.ncbins;
end
if isempty(file1), error(sprintf('Empty file %s', file1)); end
if isempty(file2), error(sprintf('Empty file %s', file2)); end
if ~exist(file1), error(sprintf('Could not find file %s', file1)); end
if ~exist(file2), error(sprintf('Could not find file %s', file2)); end

dataFile1 = loadtxt(file1);
dataFile2 = loadtxt(file2, 'delim', 9);

% decode bin file
% ---------------
bininfo = [];
metabin = [];
bCount = 1;
mCount = 1;
for iLine = 1:size(dataFile2,1)
    begPar = find(dataFile2{iLine} == '(');
    endPar = find(dataFile2{iLine} == ')');
    indEq  = find(dataFile2{iLine} == '=');
    labelInd = strfind( dataFile2{iLine}, 'label');
    if length(indEq) < 1
        error(sprintf('Could not find equal for line %d in file %s', iLine, file2));
    end
    label    = dataFile2{iLine}(labelInd+5:end);
    bin      = recodebin(dataFile2{iLine}(1:indEq(1)-1));
    if isempty(labelInd)
        error(sprintf('Could not find "label" for line %d in file %s', iLine, file2));
    end

    if length(begPar) < 1 || length(endPar) < 1 || begPar(1) > labelInd
        binrange = deblank(dataFile2{iLine}(indEq+1:labelInd-1));
        indMinus = find(binrange == '-');
        if length(indEq) < 1
            error(sprintf('Could not find minus in bin range for line %d in file %s', iLine, file2));
        end
        bin1 = recodebin(binrange(1:indMinus-1));
        bin2 = recodebin(binrange(indMinus+1:end));
        metabin(mCount).bin   = bin;
        metabin(mCount).binnum = str2num(bin(4:end));
        metabin(mCount).label = label;
        metabin(mCount).binrange = { bin1 bin2 };
        mCount = mCount+1;
    else
        bininfo(bCount).bin    = bin;
        bininfo(bCount).binnum = str2num(bin(4:end));
        bininfo(bCount).label  = label;
        bininfo(bCount).eventrange = str2num(dataFile2{iLine}(begPar(1)+1:endPar(1)-1));
        bCount = bCount+1;
    end
end

% scan events
% -----------
oldEvents = [dataFile1{:,1}]; % event type numerical
eventInds = [dataFile1{:,3}];
newName   = dataFile1(:,2);
for iEvent = 1:length(EEG.event)
    type = EEG.event(iEvent).type;
    typeInd = find(type == oldEvents);
    if length(typeInd) ~= 1
        error('Could not find type in event list');
    end
    typeIndForBin = eventInds(typeInd);
    EEG.event(iEvent).newtype  = recodeeventlabel(newName{typeInd});
    EEG.event(iEvent).eventind = typeIndForBin;
end
if ncbins
    for iEvent = 1:length(EEG.event)
        typeIndForBin = EEG.event(iEvent).eventind;
        for iBin = 1:length(bininfo)
            if any(typeIndForBin == bininfo(iBin).eventrange)
                 EEG.event(iEvent).(recodelabel(bininfo(iBin).label)) = 1;
            else EEG.event(iEvent).(recodelabel(bininfo(iBin).label)) = 0;
            end
        end
    end
end

% scan metabin and transform to bin
% ---------------------------------
wb = warning('backtrace');
warning('backtrace', 'off');
for iMeta = 1:length(metabin)
    allBins = [ bininfo.binnum ];
    bin1 = str2num(metabin(iMeta).binrange{1}(4:end));
    bin2 = str2num(metabin(iMeta).binrange{2}(4:end));
    indBin1 = find(allBins == bin1);
    indBin2 = find(allBins == bin2);
    if length(indBin1) ~= 1 || length(indBin2) ~= 1
        warning([ 'Cannot calculate contrast for ' metabin(iMeta).bin ' using ' metabin(iMeta).binrange{1} ' which is already a contrast bin' ]); 
    else
        if length(indBin2) ~= 1, error([ metabin(iMeta).binrange{2} ' not found or found twice' ]); end
        if ~isempty(intersect(bininfo(indBin1).eventrange, bininfo(indBin2).eventrange))
            warning([ 'cannot calculate contrast for ' metabin(iMeta).bin ' because ' metabin(iMeta).binrange{1} ' and ' metabin(iMeta).binrange{2} ' share some common event codes' ]);
        else
            for iEvent = 1:length(EEG.event)
                type = EEG.event(iEvent).eventind;
                if any(type == bininfo(indBin1).eventrange)
                    EEG.event(iEvent).(recodelabel(metabin(iMeta).label)) = doubledeblank(bininfo(indBin1).label);
                elseif any(type == bininfo(indBin2).eventrange)
                    EEG.event(iEvent).(recodelabel(metabin(iMeta).label)) = doubledeblank(bininfo(indBin2).label);
                else
                    EEG.event(iEvent).(recodelabel(metabin(iMeta).label)) = 'none of these';
                end
            end
        end
    end
end
warning(wb);

% return the string command
% -------------------------
com = sprintf('EEG = pop_importerplab( EEG, ''%s'',''%s'', %d);', file1, file2, ncbins);

function str = recodebin(str)

str = lower(str);
str = doubledeblank(str);
if str(1) ~= 'b'
    str = [ 'bin' str ];
elseif str(2) ~= 'i'
    str = [ 'bin' str(2:end) ];
end

function str = recodeeventlabel(str)
str(find(str == '"')) = [];

function str = recodelabel(str)

str = doubledeblank(str);
str(find(str == ' ')) = '_';
str(find(str == ')')) = '_';
str(find(str == '(')) = '_';
str(find(str == '/')) = '_';
str(find(str == '-')) = '_';
str(find(str == '&')) = '_';
if ~isempty(str2num(str(1))) str = [ 'f' str ]; end

function str = doubledeblank(str)

str = deblank(str(end:-1:1));
str = deblank(str(end:-1:1));
