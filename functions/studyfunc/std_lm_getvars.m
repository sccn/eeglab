% std_lm_getvars() - Retrieve categorical or continuous variables from a
%                    design in the STUDY structure to build the regressors
%
% Usage:
%   >>  [var_matrix,catvar_info] =
%       std_lm_getvars(STUDY,'S01','design_indx',1);
%   >>  [var_matrix,catvar_info] =
%       std_lm_getvars(STUDY,'S01','design_indx',1,'vartype','cat');
%
% Inputs:
%      STUDY       - studyset structure containing some or all files in ALLEEG
%      ALLEEG      - vector of loaded EEG datasets
%      subject     - String with the subject identifier
%
% Optional inputs:
%      vartype      - Categoricals ('cat') or continuous ('cont')
%      design_indx  - Index of the design in he STUDY structure
%
% Outputs:
% var_matrix   - Variables retrieved from the design specified in the STUDY.
%                Each column represents a factor and each row the index of
%                the variables in STUDY.design.variable.value
%                By default NaNs values will be inserted if no info from
%                that variable is not found in the original ALLEEG index,
%                This is to keep the original number of trials from ALLEEG
% catvar_info  - Information of the trials that contains all the events requested.
%                If one of the events requested is not contained in the trial,
%                then this info (second argument) would not be returned for
%                that especific trial. Notice that might be possible with
%                categotical variables, that at least one event is not
%                contained in every trial (for all of them), then this output
%                will be returned as empty.
% See also:
%
% Author: Ramon Martinez-Cancino, SCCN, 2015
%
% Copyright (C) 2015  Ramon Martinez-Cancino,INC, SCCN
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

function [var_matrix,catvar_info] = std_lm_getvars(STUDY,subject,varargin)

% Prevent empty output
var_matrix  = [];
catvar_info = [];

%% Varargin stuff
%  --------------
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end
catch
    error('std_lm_getcatvars() error: calling convention {''key'', value, ... } error'); return;
end

try, g.design;         catch, g.design      = 1 ;       end; % By default will use the first design if not especified
try, g.vartype;        catch, g.vartype     = 'cat';    end; % 'cat' or 'cont'

%% cat/cont defs
%  -------------
if strcmp(g.vartype,'cat')
    vartype = 'categorical';
elseif strcmp(g.vartype,'cont')
    vartype = 'continuous';
end
%% Checking if design
%  ------------------
if g.design > size(STUDY.design,2)
    error('std_lm_getvars() error: Invalid design index');
end

%% Checking setindex and valid subject
%  -----------------------------------
CurrentSubIndxDataset   = find(strcmp({STUDY.datasetinfo.subject},subject));
if isempty(CurrentSubIndxDataset)
    error('std_lm_getcatvars() error: A valid subject must be provided');
end
g.setindx = CurrentSubIndxDataset;

%% Getting factors
%  ---------------
varindx   = find(strcmp({STUDY.design(g.design).variable.vartype},vartype));
for i = 1:length(varindx)
    if strcmp(g.vartype,'cat')
        if ~iscell(STUDY.design(g.design).variable(varindx(i)).value{1})
            g.factors{i} = STUDY.design(g.design).variable(varindx(i)).value;
        else
            g.factors{i} = STUDY.design(g.design).variable(varindx(i)).value{:};
        end
    else
        g.factors{i} = STUDY.design(g.design).variable(varindx(i)).label;
    end
end
% Cleaning  'g.factors' from empty cells
for i = 1 : length(g.factors)
    if isempty(g.factors{i}) || all(strcmp(g.factors{i},'')), g.factors(i) = []; end
end

%% Number of trials and index
%  --------------------------
NbTrials   = 0;
dsetvect = [];
indstrialscont = [];
for i = 1 : length(g.setindx)
    nbtrials_tpm = length(STUDY.datasetinfo(g.setindx(i)).trialinfo);
    NbTrials = NbTrials + nbtrials_tpm;
    % ---
    if i == 1,
        StartEndIndx{i} = 1: nbtrials_tpm;
    else
        StartEndIndx{i} = StartEndIndx{i-1}(end) + 1 : StartEndIndx{i-1}(end) + nbtrials_tpm;
    end
    dsetvect       = [dsetvect g.setindx(i)*ones(1,length(StartEndIndx{i}))];
    indstrialscont = [indstrialscont 1: nbtrials_tpm];
end

%% Getting categorical/continuous variables from STUDY design
%  ----------------------------------------------------------
var_matrix = nan(NbTrials,length(g.factors));                   % Initializing Categorical Variables

%  Retrieving all trials and values for this subject
trialinfo  = std_combtrialinfo(STUDY.datasetinfo, g.setindx);   % Combining trialinfo
ntrials = 0;
for i = 1 : length(g.setindx)
    startendindx(i,1) = ntrials + 1;
    ntrials = ntrials + length(STUDY.datasetinfo(g.setindx(i)).trialinfo);
    startendindx(i,2) = ntrials;
end

%  Loop per variable
for i = 1 : length(varindx)
    
    % case for continuous variables
    if strcmp(g.vartype, 'cont')
        varlength = 1;
        catflag = 0;
    else
        varlength = length( STUDY.design(g.design).variable(varindx(i)).value);
        catflag = 1;
    end
    
    % Loop per Variable values
    for j = 1 : varlength
        if catflag
            if isnumeric(STUDY.design(g.design).variable(varindx(i)).value{j})
                facval      = cell2mat(STUDY.design(g.design).variable(varindx(i)).value(j));
                if length(facval)==1
                    facval_indx = find(facval == cell2mat(STUDY.design(g.design).variable(varindx(i)).value));
                else
                    facval_indx = j;
                end
            else
                if ~iscell(STUDY.design(g.design).variable(varindx(i)).value{1})
                    facval = cell2mat(STUDY.design(g.design).variable(varindx(i)).value(j));
                    facval_indx = find(strcmp(facval,STUDY.design(g.design).variable(varindx(i)).value));
                else
                    facval_indx = j;
                end
                
            end
        end
        
        % No loop per dataset since we merged datasetinfo
        if catflag
            if isnumeric(STUDY.design(g.design).variable(varindx(i)).value{j})
                varval = cell2mat(STUDY.design(g.design).variable(varindx(i)).value(j));
            else
                if ~iscell(STUDY.design(g.design).variable(varindx(i)).value{1})
                    varval = STUDY.design(g.design).variable(varindx(i)).value(j);
                else
                    varval = STUDY.design(g.design).variable(varindx(i)).value{j};
                end
            end
        else
            varval = '';
        end
        [trialindsx, eventvals] = std_gettrialsind(trialinfo,STUDY.design(g.design).variable(varindx(i)).label, varval);
        if ~isempty(trialindsx)
            % case for continuous variables
            if ~catflag
                facval_indx = eventvals;
            end
            var_matrix(trialindsx,i) = facval_indx;
        end
    end
end

%% Getting trialindex  and names for overlapped conditions 
% This means for the trials that contains all the type of events requested.
% If one of the events requested is not contained in the trial, then this info (second argument) would not be returned
%  -------------------------------------------------------
tmpmat = var_matrix;
tmpindx = find(isnan(tmpmat));
[I,tmp] = ind2sub(size(tmpmat),tmpindx); clear tmp; %#ok<ASGLU>
tmpmat(I,:) = [];
if ~isempty(tmpmat)
%     If all the rows in 'var_matrix' contains at least one NaN, then this
%     loop will no be executed
    if strcmp(g.vartype,'cat')
        ind = 1;
        comb = unique(tmpmat,'rows');
        for i = 1: size(comb,1)
            TrialIndx_tmp{i} = find(sum(repmat(comb(i,:),[size(var_matrix,1),1]) == var_matrix,2) == size(var_matrix,2));
            tmpsets = dsetvect(TrialIndx_tmp{i});
            uniqtmpset = unique(tmpsets);
            if length(uniqtmpset) ~= 1
                tmpvect = dsetvect;
                tmpvect(setdiff(1:length(dsetvect),TrialIndx_tmp{i})) = 0;
                
                for k = 1: length(uniqtmpset)
                    TrialIndx_datasetinfo{ind} = find(uniqtmpset(k) == tmpvect);
                    datasets{ind} = uniqtmpset(k);
                    ind = ind + 1;
                end
            else
                TrialIndx_datasetinfo{ind} = TrialIndx_tmp{i};
                datasets{ind}              = uniqtmpset;
            end
            ind = ind + 1;
        end
    else
        datasets = g.setindx;
        TrialIndx_datasetinfo = StartEndIndx;
        if ~isempty(tmpindx)
            for i = 1:length(TrialIndx_datasetinfo)
                for j = 1:length(I)
                    if ismember(I(j),TrialIndx_datasetinfo{i}), TrialIndx_datasetinfo{i}(I(j)) = []; end
                end
            end
        end
    end
    
    for i = 1:length(TrialIndx_datasetinfo)
        TrialIndx_sets{i} =  indstrialscont(TrialIndx_datasetinfo{i});
    end
    
    %% Outputs
    %  -------
    catvar_info.datasetinfo_trialindx   = TrialIndx_sets;  % Indices at datasetinfo.trialinfo {1:Ntrials1} {1:Ntrials2}
    catvar_info.concat_trialindx        = StartEndIndx;    % Indices at [ 1 : Ntrials1 , Ntrials1+1 : Ntrials2]
    catvar_info.datasetinfo_concatindx  = g.setindx;
    catvar_info.dataset                 = datasets;
end
