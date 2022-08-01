% pop_fusechanrej() - Make sure the same subject and session have the same
%                     removed. If not remove channel not in common. 
% Usage:
%   >>  ALLEEG = pop_fusechanrej(ALLEEG);
%
% see also: pop_clean_rawdata

% Copyright (C) 2022 Arnaud Delorme, UCSD
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

function EEG = pop_fusechanrej(EEG, varargin)

if nargin < 1
    help pop_fusechanrej;
    return;
end

if length(EEG) ==  1
    return;
end

% find common datasets
% same code below as for pop_runica
% ---------------------------------
allsubjects = { EEG.subject };
allsessions = { EEG.session };
alltags     = zeros(1,length(allsubjects));
if any(cellfun('isempty', allsubjects))
    warning( [ 'Cannot fuse channel rejection within subject because' ...
        'subject names missing from at least one dataset file.' 10 ...
        'Subject names must be stored within the datasets. To do so,' 10 ...
        'use the STUDY > Edit STUDY Info menu and check the box' 10 ...
        '"Dataset info (condition, group, ...) differs from study info..."' ]);
    return;
end
dats = {};
for index = 1:length(allsubjects)
    if ~alltags(index)
        allinds = strmatch(allsubjects{index}, allsubjects, 'exact');
        rmind = [];
        % if we have different sessions they will not be concatenated
        for tmpi = setdiff_bc(allinds,index)'
            if ~isequal(allsessions(index), allsessions(tmpi))
                rmind = [rmind tmpi];
            end
        end
        allinds = setdiff_bc(allinds, rmind);
        fprintf('Found %d datasets for subject ''%s'' session %d\n', length(allinds), allsubjects{index}, allsessions{index});
        dats = { dats{:} allinds };
        alltags(allinds) = 1;
    end
end

eeglab_options;
for index = 1:length(dats)
    if length(dats{index}) == 1
        TMPALLEEG = EEG(dats{index}(1));
        if ~isempty(TMPALLEEG(1).session)
            fprintf('Skipping selecting common channels accross datasets for subject %s session %s (only 1 dataset)\n', TMPALLEEG(1).subject, num2str(TMPALLEEG(1).session));
        else
            fprintf('Skipping selecting common channels accross datasets for subject %s (only 1 dataset)\n', TMPALLEEG(1).subject);
        end
    else
        TMPALLEEG = EEG(dats{index});
        commonChans = myintersect(TMPALLEEG.chanlocs);
        if ~isempty(TMPALLEEG(1).session)
            fprintf('Selecting common channels accross datasets for subject %s session %s ***************\n', TMPALLEEG(1).subject, num2str(TMPALLEEG(1).session));
        else
            fprintf('Selecting common channels accross datasets for subject %s ***************\n', TMPALLEEG(1).subject);
        end
        for iSet = 1:length(TMPALLEEG)
            if length(commonChans) ~= TMPALLEEG(iSet).nbchan
                TMPEEG = eeg_retrieve(TMPALLEEG, iSet);
                TMPEEG = pop_select(TMPEEG, 'channel', { commonChans.labels });
                TMPEEG = eeg_checkset(TMPEEG);
                TMPEEG.saved = 'no';
                if option_storedisk
                    TMPEEG = pop_saveset(TMPEEG, 'savemode', 'resave');
                    TMPEEG = update_datafield(TMPEEG);
                end
                EEG = eeg_store(EEG, TMPEEG, dats{index}(iSet));
                if option_storedisk
                    EEG(dats{index}(iSet)).saved = 'yes'; % eeg_store by default set it to no
                end
            end
        end
    end
end

% same as in eeg_eval
% -------------------
function EEG = update_datafield(EEG)
    if ~isfield(EEG, 'datfile'), EEG.datfile = ''; end
    if ~isempty(EEG.datfile)
        EEG.data = EEG.datfile;
    else 
        EEG.data = 'in set file';
    end
    EEG.icaact = [];    

% without losing the order information
% ---------------------------------------
function alllocs = myintersect(locs1, locs2, varargin)

    if length(varargin) >= 1
        tmplocs = myintersect(locs1, locs2);
        alllocs = myintersect(tmplocs, varargin{1}, varargin{2:end});
        return
    end

   labs1 = { locs1.labels };
   labs2 = { locs2.labels };
   
   count1 = 1;
   count2 = 1;
   count3 = 1;
   alllocs = locs1; alllocs(:) = [];
   while count1 <= length(locs1) && count2 <= length(locs2)
       
       if strcmpi(labs1{count1}, labs2{count2}) % same label
           alllocs(count3) = locs1(count1); % copy
           count1 = count1 + 1;
           count2 = count2 + 1;
           count3 = count3 + 1;
       elseif isempty(strmatch(labs1{count1}, labs2, 'exact'))
           count1 = count1 + 1;
       elseif isempty(strmatch(labs2{count2}, labs1, 'exact'))
           count2 = count2 + 1;
       else 
           count1 = count1 + 1;
           count2 = count2 + 1;
       end
   end
   return   

   a(1).labels = '1';
   a(2).labels = '3';
   a(3).labels = '4';
   a(4).labels = '7';
   
   b(1).labels = '3';
   b(2).labels = '4';
   b(3).labels = '5';
   b(4).labels = '6';

   c(1).labels = '1';
   c(2).labels = '2';
   c(3).labels = '4';
   c(4).labels = '7';
   
   myintersect_delete(a,b)
   myintersect_delete(a,b,c)