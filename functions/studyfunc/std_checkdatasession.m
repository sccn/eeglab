% std_checkdatasession() - Check/fill session field in STUDY. Based in the
% name and IC decomposition of datases, the function determine which
% ones are coming from the same sessions. 
% It will assign the same index number to the sets from the same subject coming
% from the same session.
%
% Usage:
%   >>  clustinfo = std_checkdatasession(STUDY, ALLEEG);
%
% Inputs:
%      STUDY
%      ALLEEG   - vector of loaded EEG datasets
%
% Optional inputs:
%
% session          - vector of the same size of STUDY.datasetinfo with with
%                    session index. Default empty
% verbose          - [0,1] Default 1
%
%
% Outputs:
%    STUDY  - studyset structure containing some or all files in ALLEEG
%             (datasetinfo.session updated)
%    flags  - Vector of the dimension of the numbre of subjects (Ordered as
%    in STUDY.datasetinfo) with the information if the subject have dataset
%    from different sesssions (1) or not (0)
%
% See also:
%   std_plotinfocluster
%
% Author: Ramon Martinez-Cancino, SCCN, 2014
%
% Copyright (C) 2014  Ramon Martinez-Cancino,INC, SCCN
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

function [STUDY,flags] = std_checkdatasession(STUDY, ALLEEG,varargin)


%--------------------------------------------------------------------------
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end
catch
    disp('std_checkdatasession() error: calling convention {''key'', value, ... } error'); return;
end

try, g.session;    catch, g.session   = [];        end; % Index for sessions
try, g.verbose;    catch, g.verbose    = 1;        end;   % By default verbose

isession = ~(cellfun(@isempty, {STUDY.datasetinfo.session}));
if std(isession) == 0
    flags = zeros(size(isession));
    return; 
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
UniqueSubj = unique({STUDY.datasetinfo.subject});

if isempty(g.session)
    
    % Step 1: checking IC decomposition
    alleeg_indx  = cell2mat({STUDY.datasetinfo.index});         % data index from study
    sameica = std_findsameica(ALLEEG);
    
    if length(sameica) ~= length(alleeg_indx)
        if g.verbose, display('--- Data from same sessions found in the STUDY ---'); end
        
        UniqueSubj = unique({STUDY.datasetinfo.subject});
        
        % Step 2: Identify same ica decomposition
        for i = 1 : length(UniqueSubj)
            SubjInd = find(ismember({STUDY.datasetinfo.subject},UniqueSubj(i)));
            temp_sameica = std_findsameica(ALLEEG(SubjInd));      % Finding same ica
            
            c = 1;
            for j = 1:length(temp_sameica)
                [STUDY.datasetinfo(SubjInd(cell2mat(temp_sameica(j)))).session] = deal(c);
                if g.verbose, display(['Datasets with indices [' num2str(SubjInd(cell2mat(temp_sameica(j)))) '] assigned to session ' num2str(c)]); end
                c = c + 1;
                
            end
        end
        
        
    else
        [STUDY.datasetinfo.session] = deal(1);
    end
    
    if g.verbose, display('--- Session fields in current STUDY succesfully updated ---'); end
    
    flags = zeros(1,length(UniqueSubj));
    for i = 1: length(UniqueSubj)
        if length({STUDY.datasetinfo(find(strcmp({STUDY.datasetinfo.subject},UniqueSubj{i}))).session}) ~= 1
            flags(i) = ~isequal(STUDY.datasetinfo(find(strcmp({STUDY.datasetinfo.subject},UniqueSubj{i}))).session);
        end
    end
    
% Updating field with custom info
elseif ~(isempty(g.session))
    
    if length(g.session(:))== length(STUDY.datasetinfo)
        
        % I dont like this loop... but it will be here until I find
        % something more elegant
        for i = 1:length(STUDY.datasetinfo)
            STUDY.datasetinfo(i).session = g.session(i); %Assigning values
        end
        
        if g.verbose, display('--- Session fields in current STUDY succesfully updated ---'); end
        
    else
        error('Error in std_checkdatasession(): Vector dimensions inconsisten with STUDY');
        eeglab_error;
        return
    end
end
