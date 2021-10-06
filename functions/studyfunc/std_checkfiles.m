% std_checkfiles() - Check all STUDY files consistency
%
% Usage:    
%                >> boolval = std_checkfiles(STUDY, ALLEEG);
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - All EEGLAB datasets
%
% Outputs:
%   boolval    - [0|1] 1 if uniform
%
% Authors: Arnaud Delorme, CERCO, 2010-

% Copyright (C) Arnaud Delorme, CERCO
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

function [boolval npersubj] = std_checkfiles(STUDY, ALLEEG);

if nargin < 2
    help std_checkfiles;
    return;
end
return;

filetypes = { 'daterp' 'datspec' 'datersp' 'datitc' 'dattimef' ...
              'icaerp' 'icaspec' 'icaersp' 'icaitc' 'icatopo' };

% set channel interpolation mode
% ------------------------------
uniformchannels = std_uniformsetinds( STUDY );

disp('---------------------------------------------');
disp('Checking data files integrity and consistency');
passall = 1;
for index = 1:length(filetypes)
    
    % scan datasets
    % -------------
    [ tmpstruct compinds filepresent ] = std_fileinfo(ALLEEG, filetypes{index});
    if ~isempty(tmpstruct)
        fprintf('Files of type "%s" detected, checking...',  filetypes{index});
    end
    
    % check if the structures are equal
    % ---------------------------------
    notequal = any(~filepresent);
    if ~isempty(tmpstruct)
        if any(filepresent == 0)
            fprintf(' Error, some files inconsistent or missing\n');
            notequal = 1;
            passall = 0;
        else
            firstpass = 1;
            fields = fieldnames(tmpstruct);
            for f_ind = 1:length(fields)
                
                firstval = getfield(tmpstruct, {1}, fields{f_ind});
                for dat = 2:length(tmpstruct)
                    tmpval = getfield(tmpstruct, {dat}, fields{f_ind});
                    if ~isequal(firstval, tmpval)
                        
                        % check for NaNs
                        if iscell(firstval) && iscell(tmpval)
                            for cind = 1:length(firstval)
                                if isreal(firstval{cind}) && ~isempty(firstval{cind}) && isnan(firstval{cind}(1))
                                    firstval{cind} = 'NaN';
                                end
                            end
                            for cind = 1:length(tmpval)
                                if isreal(tmpval{cind}) && ~isempty(tmpval{cind}) && isnan(tmpval{cind}(1))
                                    tmpval{cind} = 'NaN';
                                end
                            end
                        end
                        
                        if ~isequal(firstval, tmpval)
                            if ~strcmpi(fields{f_ind}, 'labels') || strcmpi(uniformchannels, 'on')
                                if firstpass == 1, fprintf('\n'); firstpass = 0; end
                                fprintf('  Error, difference across data files for field "%s"\n', fields{f_ind});
                                notequal = 1;
                                passall = 0;
                                break;
                            end
                        end
                    end
                end
            end
        end
    end
    
    % check the consistency of changrp and cluster with saved information
    % -------------------------------------------------------------------
    if isempty(tmpstruct), notequal = 1; end
    if filetypes{index}(1) == 'd' && notequal == 0        
        % scan all channel labels
        % -----------------------
        if isfield(tmpstruct(1), 'labels')
            for cind = 1:length(STUDY.changrp)
                if notequal == 0
                    for inddat = 1:length(ALLEEG)
                        tmpind = cellfun(@(x)(find(x == inddat)), STUDY.changrp(cind).setinds(:), 'uniformoutput', false);
                        indnonempty = find(~cellfun(@isempty, tmpind(:)));
                        if ~isempty(indnonempty)
                            tmpchan = STUDY.changrp(cind).allinds{indnonempty}(tmpind{indnonempty}); % channel index for dataset inddat

                            tmpchan2 = strmatch(STUDY.changrp(cind).name, tmpstruct(inddat).labels, 'exact'); % channel index in file
                            if ~isempty(tmpchan2) || ~strcmpi(filetypes{index}, 'datspec') % the last statement is because channel labels did not use to be saved in spec files
                                if ~isequal(tmpchan2, tmpchan)
                                    fprintf('\nError: channel index in STUDY.changrp(%d) for dataset %d is "%d" but "%d" in data files\n', cind, inddat, tmpchan, tmpchan2);
                                    notequal = 1;
                                    break;
                                end
                            end
                        end
                    end
                end
            end
        end
    elseif notequal == 0 && ~isempty(STUDY.cluster) % components
        % check that the cell structures are present
        % ------------------------------------------
        if ~isfield(STUDY.cluster, 'setinds')
            STUDY.cluster(1).setinds = [];
            STUDY.cluster(1).allinds = [];
        end
        for cind = 1:length(STUDY.cluster)
            if isempty(STUDY.cluster(cind).setinds)
                STUDY.cluster(cind) = std_setcomps2cell(STUDY, cind);
            end
        end
        for cind = 1:length(STUDY.cluster)
            if notequal == 0
                for inddat = 1:length(ALLEEG)
                    tmpind = cellfun(@(x)(find(x == inddat)), STUDY.cluster(cind).setinds(:), 'uniformoutput', false);
                    indnonempty = find(~cellfun(@isempty, tmpind(:)));
                    tmpcomp = [];
                    for jind = 1:length(indnonempty)
                        tmpcomp = [ tmpcomp STUDY.cluster(cind).allinds{indnonempty(jind)}(tmpind{indnonempty(jind)}) ];
                    end
                    
                    if ~isempty(setdiff(tmpcomp, compinds{inddat}))
                        if ~(isempty(compinds{inddat}) && strcmpi(filetypes{index}, 'icatopo'))
                            fprintf('\nError: some components in clusters %d are absent from .%s files\n', cind, filetypes{index});
                            notequal = 1;
                            passall  = 0;
                            break;
                        end
                    end
                end
            end
        end
    end
    if notequal == 0, fprintf(' Pass\n'); end
end

if ~passall
    disp('**** Recompute any measure above not receiving a "Pass" by')
    disp('**** calling menu items "STUDY > Precompute Channel/Component measures" ');
    disp('**** and by selecting the "Recompute even if present on disk" checkbox');
end
disp('Checking completed.');
