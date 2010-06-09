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

function [ factor factorvals subjects ] = std_getindvar(STUDY, mode, scandesign)

if nargin < 1
    help std_getindvar;
    return;
end;

factor     = {};
factorvals = {};
subjects   = {};
if nargin < 2, mode = 'both'; end;
if nargin < 3, scandesign = 0; end;
countfact = 1;
setinfo = STUDY.datasetinfo;

if strcmpi(mode, 'datinfo') || strcmpi(mode, 'both')
    % get trial info
    ff = fieldnames(setinfo);
    ff = setdiff(ff, { 'filepath' 'filename' 'subject' 'index' 'ncomps' 'comps' 'trialinfo' });
    for index = 1:length(ff)
        if isstr(getfield(setinfo(1), ff{index}))
            eval( [ 'tmpvals = unique({ setinfo.' ff{index} '});' ] );
            if length(tmpvals) > 1
                factor{    countfact} = ff{index};
                factorvals{countfact} = tmpvals;
                
                % get subject for each factor value
                for c = 1:length(tmpvals)
                    eval( [ 'datind = strmatch(tmpvals{c}, { setinfo.' ff{index} '}, ''exact'');' ] );
                    subjects{  countfact}{c} = unique( { setinfo(datind).subject } );
                end;
                countfact = countfact + 1;
            end;
        end;
    end;
end;

% add ind. variables for trials
% -----------------------------
if strcmpi(mode, 'trialinfo') || strcmpi(mode, 'both')
    % add trial info
    if isfield(setinfo, 'trialinfo')
        ff = fieldnames(setinfo(1).trialinfo);
        for index = 1:length(ff)
            if isstr(getfield(setinfo(1).trialinfo(1), ff{index}))
                alltmpvals = {};
                for ind = 1:length(setinfo)
                    if isfield(setinfo(ind).trialinfo, ff{index})
                         eval( [ 'tmpvals = unique({ setinfo(ind).trialinfo.' ff{index} ' });' ] );
                    else tmpvals = {};
                    end;
                    alltmpvals = { alltmpvals{:} tmpvals{:} };
                end;
                alltmpvals = unique(alltmpvals);
                if length(alltmpvals) > 1
                    factor{    countfact} = ff{index};
                    factorvals{countfact} = alltmpvals;
                    subjects{  countfact} = {};
                    countfact = countfact + 1;
                end;
            end;
        end;
    end;
end;

% scan existing design for additional combinations
% ------------------------------------------------
if scandesign
    for desind = 1:length(STUDY.design)
        pos1 = strmatch(STUDY.design(desind).variable(1).label, factor, 'exact');
        pos2 = strmatch(STUDY.design(desind).variable(2).label, factor, 'exact');
        if ~isempty(pos1), add1 = setdiff(STUDY.design(desind).variable(1).value, factorvals{pos1}); else add1 = []; end;
        if ~isempty(pos2), add2 = setdiff(STUDY.design(desind).variable(2).value, factorvals{pos2}); else add2 = []; end;
        if ~isempty(add1), factorvals{pos1} = { factorvals{pos1}{:} add1{:} }; end;
        if ~isempty(add2), factorvals{pos2} = { factorvals{pos2}{:} add2{:} }; end;
    end;
end;
