% std_getindvar - get independent variables of a STUDY
%
% Usage:
%   [indvar indvarvals] = std_getindvar(STUDY);
%   [indvar indvarvals] = std_getindvar(STUDY, mode, scandesign);
%
% Input:
%   STUDY - EEGLAB STUDY structure
%   mode  - ['datinfo'|'trialinfo'|'both'] get independent variables
%           linked to STUDY.datasetinfo, STUDY.datasetinfo.trialinfo or
%           both. Default is 'both'.
%   scandesign - [0|1] scan STUDY design for additional combinations of
%           independent variable values. Default is 0.
%
% Output:
%   indvar     - [cell array] cell array of independent variable names 
%   indvarvals - [cell array] cell array of independent variable values
%
% Authors: A. Delorme, CERCO/CNRS and SCCN/UCSD, 2010-

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2010, arno@sccn.ucsd.edu
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

function [ factor, factorvals, subjects, paired ] = std_getindvar(STUDY, mode, scandesign)

if nargin < 1
    help std_getindvar;
    return;
end

factor     = {};
factorvals = {};
paired     = {};
subjects   = {};
if nargin < 2, mode = 'both'; end
if nargin < 3, scandesign = 0; end
countfact = 1;
setinfo = STUDY.datasetinfo;

if strcmpi(mode, 'datinfo') || strcmpi(mode, 'both')
    % get trial info
    ff = fieldnames(setinfo);
    ff = setdiff_bc(ff, { 'filepath' 'filename' 'subject' 'index' 'ncomps' 'comps' 'trialinfo' });
    for index = 1:length(ff)
        % check we have the same data type
        allSetinfoVals = { setinfo.(ff{index}) };
        if any(cellfun(@ischar, allSetinfoVals))
            if ~all(cellfun(@ischar, allSetinfoVals))
                allSetinfoVals = cellfun(@num2str, allSetinfoVals, 'uniformoutput', false);
            end
        end

        if ischar(allSetinfoVals{1})
            tmpvals = unique_bc(allSetinfoVals);
            if length(tmpvals) > 1
                factor{    countfact} = ff{index};
                factorvals{countfact} = tmpvals;

                % get subject for each factor value
                intersectSubject = { setinfo(:).subject };
                for c = 1:length(tmpvals)
                    eval( [ 'datind = strmatch(tmpvals{c}, { setinfo.' ff{index} '}, ''exact'');' ] );
                    subjects{  countfact}{c} = unique_bc( { setinfo(datind).subject } );
                    intersectSubject = intersect(intersectSubject, subjects{  countfact}{c});
                end
                
                numValues = cellfun(@length, subjects{  countfact});
                if length(intersectSubject) == numValues(1) && length(unique(numValues)) == 1
                    paired{countfact} = 'on'; % trial data always paired
                else
                    paired{countfact} = 'off';
                end
                countfact = countfact + 1;
            end
        else
            tmpvals = unique_bc([ allSetinfoVals{:} ] );
            if length(tmpvals) > 1
                factor{    countfact} = ff{index};
                factorvals{countfact} = mattocell(tmpvals);

                % get subject for each factor value
                intersectSubject = { setinfo(:).subject };
                for c = 1:length(tmpvals)
                    eval( [ 'datind = find(tmpvals(c) == [ setinfo.' ff{index} ']);' ] );
                    subjects{  countfact}{c} = unique_bc( { setinfo(datind).subject } );
                    intersectSubject = intersect(intersectSubject, subjects{  countfact}{c});
                end
                
                numValues = cellfun(@length, subjects{  countfact});
                if length(intersectSubject) == numValues(1) && length(unique(numValues)) == 1
                    paired{countfact} = 'on'; % trial data always paired
                else
                    paired{countfact} = 'off';
                end
                countfact = countfact + 1;
            end
        end
    end
end

% add ind. variables for trials
% -----------------------------
if strcmpi(mode, 'trialinfo') || strcmpi(mode, 'both')
    % add trial info
    if isfield(setinfo, 'trialinfo')
        ff = fieldnames(setinfo(1).trialinfo);
        for index = 1:length(ff)
            % check if any of the datasets are using string for event type
            allFieldsPresent = cellfun(@(x)(isfield(x, ff{index})), { setinfo.trialinfo });                
            allFirstVal = cellfun(@(x)(getfield(x, ff{index})), { setinfo(allFieldsPresent).trialinfo }, 'uniformoutput', false);

            if any(cellfun(@isstr, allFirstVal))
                alltmpvals = {};
                for ind = 1:length(setinfo)
                    if isfield(setinfo(ind).trialinfo, ff{index})
                        eval( [ 'tmpTrialVals = { setinfo(ind).trialinfo.' ff{index} ' };' ] );
                        if isnumeric(tmpTrialVals{1})
                            % convert to string if necessary
                            tmpTrialVals = cellfun(@num2str, tmpTrialVals, 'uniformoutput', false);
                        end
                        tmpvals = unique_bc(tmpTrialVals);
                    else tmpvals = {};
                    end
                    if isempty(alltmpvals)
                        alltmpvals = tmpvals;
                    else
                        alltmpvals = { alltmpvals{:} tmpvals{:} };
                    end
                end
                alltmpvals = unique_bc(alltmpvals);
                if length(alltmpvals) > 1
                    factor{    countfact} = ff{index};
                    factorvals{countfact} = alltmpvals;
                    subjects{  countfact} = {};
                    paired{countfact}     = 'on'; % trial data always paired
                    countfact = countfact + 1;
                end
            else
                alltmpvals = [];
                for ind = 1:length(setinfo)
                    if isfield(setinfo(ind).trialinfo, ff{index}) && ~iscell( setinfo(ind).trialinfo(1).(ff{index}))
                        eval( [ 'tmpvals = unique_bc([ setinfo(ind).trialinfo.' ff{index} ' ]);' ] );
                    else tmpvals = [];
                    end
                    alltmpvals = [ alltmpvals tmpvals ];
                end
                alltmpvals = unique_bc(alltmpvals);
                if length(alltmpvals) > 1
                    factor{    countfact} = ff{index};
                    factorvals{countfact} = mattocell(alltmpvals);
                    subjects{  countfact} = {};
                    paired{countfact}     = 'off'; % pairing is irrelevant for continuous var
                    countfact = countfact + 1;
                end
            end
        end
    end
end

% scan existing design for additional combinations
% ------------------------------------------------
if scandesign
    for desind = 1:length(STUDY.design)
        for iVar = 1:length(STUDY.design(desind).variable)
            if strcmpi(STUDY.design(desind).variable(iVar).vartype, 'categorical')
                pos1 = strmatch(STUDY.design(desind).variable(iVar).label, factor, 'exact');
                if ~isempty(pos1), add1 = mysetdiff(STUDY.design(desind).variable(iVar).value, factorvals{pos1}); else add1 = []; end
                for addInd = 1:length(add1)
                    duplicate = 0;
                    for iVarVal = 1:length(factorvals{pos1})
                        if isequal(factorvals{pos1}{iVarVal}, add1{addInd})
                            duplicate = 1;
                        end
                    end
                    if ~duplicate
                        factorvals{pos1} = { factorvals{pos1}{:} add1{addInd} }; 
                    end
                end
            end
        end
    end
end

function cellout = mysetdiff(cell1, cell2);

    if ischar(cell2{1})
         indcell = cellfun(@iscell, cell1);
    else indcell = cellfun(@(x)(length(x)>1), cell1);
    end
    cellout = cell1(indcell);

%     if ischar(cell1{1}) && ischar(cell2{1})
%          cellout = setdiff_bc(cell1, cell2);
%     elseif ~ischar(cell1{1}) && ~ischar(cell2{1})
%         cellout = mattocell(setdiff( [ cell1{:} ], [ cell2{:} ]));
%     elseif ischar(cell1{1}) && ~ischar(cell2{1})
%          cellout = setdiff_bc(cell1, cellfun(@(x)(num2str(x)),cell2, 'uniformoutput', false));
%     elseif ~ischar(cell1{1}) && ischar(cell2{1})
%          cellout = setdiff_bc(cellfun(@(x)(num2str(x)),cell1, 'uniformoutput', false), cell2);
%     end
    
