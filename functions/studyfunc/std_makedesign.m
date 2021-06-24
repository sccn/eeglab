% std_makedesign() - create a new or edit an existing STUDY.design by 
%                    selecting specific factors to include in subsequent 
%                    1x2 or 2x2 STUDY measures and statistical computations 
%                    for this design. A STUDY may have many factors 
%                    (task or stimulus conditions, subject groups, session 
%                    numbers, trial types, etc.), but current EEGLAB
%                    STUDY statistics functions apply only to at most two 
%                    (paired or unpaired) factors. A STUDY.design may 
%                    also be (further) restricted to include only specific 
%                    subjects, datasets, or trial types.
% Usage:
%   >> [STUDY] = std_makedesign(STUDY, ALLEEG); % create a default design 
%   >> [STUDY] = std_makedesign(STUDY, ALLEEG, designind, 'key', 'val' ...);   
%
% Inputs:
%   STUDY      - EEGLAB STUDY set
%   ALLEEG     - vector of the EEG datasets included in the STUDY structure 
%   designind   - [integer > 0] index (number) of the new STUDY design {default: 1}
%
% Optional inputs:
%  'name'     - ['string'] mnemonic design name (ex: 'Targets only') 
%               {default: 'Design d', where d = designind} 
%  'variable1' - ['string'] - first independent variable or contrast. Must 
%               be a field name in STUDY.datasetinfo or in 
%               STUDY.datasetinfo.trialinfo. Typical choices include (task 
%               or other) 'condition', (subject) 'group', (subject) 'session', 
%               or other condition/group/session factors in STUDY.datasetinfo 
%               -- for example (subject factor) 'gender' or (condition factor) 
%              'timeofday', etc. If trial type variables are defined in 
%               STUDY.datasetinfo.trialinfo, they may also be used here 
%               -- for example, 'stimcolor'. However, in this case datasets 
%               consist of heterogeneous sets of trials of different types, 
%               so many dataset Plot and Tools menu items may not give 
%               interpretable results and will thus be made unavailable for 
%               selection {default: 'condition'}
%  'values1'   - {cell array of 'strings'} - 'variable1' instances to include 
%               in the design. For example, if 'variable1' is 'condition'and 
%               three values for 'condition' (e.g., 'a' , 'b', and 'c')
%               are listed in STUDY.datasetinfo, then 'indval1', { 'a' 'b' } 
%               will contrast conditions 'a' and 'b', and datasets for 
%               condition 'c' will be ignored. To combine conditions, use 
%               nested '{}'s. For example, to combine conditions 'a' and 
%               'b' into one condition and contrast it to condition 'c', 
%               specify  'indval1', { { 'a' 'b' } 'c' } {default: all values 
%               of 'variable1' in STUDY.datasetinfo}
%  'variable2' - ['string'] - second independent variable name, if any. Typically, 
%               this might refer to ('unpaired') subject group or (typically 
%               'paired') session number, etc.
%  'values2'  - {cell array of 'strings'} - variable2 values to include in the 
%               design {default: all}. Here, 'var[12]' must be field names 
%               in STUDY.datasetinfo or  STUDY.datasetinfo.trialinfo. 
%  'datselect'  - {cell array} select specific datasets and/or trials: 'datselect',
%               {'var1' {'vals'}  'var2' {'vas'}}. Selected datasets must 
%               meet all the specified conditions. For example, 'datselect',  
%               { 'condition' { 'a' 'b' } 'group' { 'g1' 'g2' } }  will 
%               select only datasets from conditions 'a' OR 'b' AND only 
%               subjects in groups 'g1' OR 'g2'. If 'subjselect' is also 
%               specified, only datasets meeting both criteria are included. 
%               'variable1' and 'variable2' will only consider
%               the values after they have passed through 'datselect' and 
%               'subjselect'. For instance, if conditions { 'a' 'b' 'c' } 
%               exist and conditions 'a' is removed by 'datselect', the only 
%               two conditions that will be considered are 'b' and 'c' 
%               (which is then equivalent to using 'variable1vals' to specify
%               values for the 'condition' factor. Calls function 
%               std_selectdataset() {default: select all datasets}
%  'subjselect' - {cell array} subject codes of specific subjects to include 
%               in the STUDY design {default: all subjects in the specified 
%               conditions, groups, etc.} If 'datselect' is also specified,
%               only datasets meeting both criteria are included. 
%  'rmfiles'  - ['on'|'off'] remove from the STUDY all data measure files 
%               NOT included in this design. Selecting this option will 
%               remove all the old measure files associated with the previous 
%               definition of this design. {default: 'off'}  
%  'filepath' - [string] file path for saving precomputed files. Default is
%               empty meaning it is in the same folder as the data.
%  'delfiles' - ['on'|'off'|'limited'] delete data files
%               associated with the design specified as parameter. 'on'
%               delete all data files related to the design. 'limited'
%               deletes all data files contained in the design. 'limited'
%               will not delete data files from other STUDY using the same
%               files. Default is 'off'.
%
% Output:
%      STUDY - The input STUDY with a new design added and designated 
%              as the current design.
%
% Examples:
%  STUDY = std_makedesign(STUDY, ALLEEG); % make default design
%
%  % create design with 1 independent variable equal to 'condition'
%  STUDY = std_makedesign(STUDY, ALLEEG, 2, 'variable1', 'condition');
%
%  % create design with 1 independent variable equal to 'condition'
%  % but only consider the sub-condition 'stim1' and 'stim2' - of course
%  % such conditions must be present in the STUDY
%  STUDY = std_makedesign(STUDY, ALLEEG, 2, 'variable1', 'condition', ...
%                    'values1', { 'stim1' 'stim2' }); 
%
%  % create design and only include subject 's1' and 's2'
%  STUDY = std_makedesign(STUDY, ALLEEG, 2, 'variable1', 'condition', ...
%                     'subjselect', { 's1' 's2' });
%
% Author: Arnaud Delorme, Institute for Neural Computation UCSD, 2010-

% Copyright (C) Arnaud Delorme, arno@salk.edu
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

function [STUDY, com] = std_makedesign(STUDY, ALLEEG, designind, varargin)

if nargin < 2
    help std_makedesign;
    return;
end
if nargin < 3 
    designind = 1;
end

defdes.name = sprintf('STUDY.design %d', designind);
defdes.cases.label = 'subject';
defdes.cases.value = {};
defdes.variable(1).label = 'condition';
defdes.variable(2).label = 'group';
defdes.variable(3).label = '';
defdes.variable(4).label = '';
defdes.variable(1).value = {};
defdes.variable(2).value = {};
defdes.variable(3).value = {};
defdes.variable(4).value = {};
defdes.variable(1).vartype = 'categorical';
defdes.variable(2).vartype = 'categorical';
defdes.filepath = '';
defdes.include = {};
orivarargin = varargin;
if ~isempty(varargin) && isstruct(varargin{1})
    defdes = varargin{1};
    varargin(1) = [];
end
if isempty(defdes.filepath), defdes.filepath = ''; end
if length(defdes.variable) == 0, defdes.variable(1).label = 'continuous'; defdes.variable(1).vartype = 'categorical'; defdes.variable(1).value = {}; end
if length(defdes.variable) == 1, defdes.variable(2).label = 'continuous'; defdes.variable(2).vartype = 'categorical'; defdes.variable(2).value = {}; end
if length(defdes.variable) == 2, defdes.variable(3).label = 'continuous'; defdes.variable(3).vartype = 'categorical'; defdes.variable(3).value = {}; end
if length(defdes.variable) == 3, defdes.variable(4).label = 'continuous'; defdes.variable(4).vartype = 'categorical'; defdes.variable(4).value = {}; end
for iVar = 1:4
    if isempty(defdes.variable(iVar).vartype), defdes.variable(iVar).vartype = 'continuous'; end
    if isempty(defdes.variable(iVar).label)  , defdes.variable(iVar).label = ''; end
end
opt = finputcheck(varargin,  {'variable1'     'string'    []     defdes.variable(1).label;
                              'variable2'     'string'    []     defdes.variable(2).label;
                              'variable3'     'string'    []     defdes.variable(3).label;
                              'variable4'     'string'    []     defdes.variable(4).label;
                              'values1'       {'real','cell' } []     defdes.variable(1).value;
                              'values2'       {'real','cell' } []     defdes.variable(2).value;
                              'values3'       {'real','cell' } []     defdes.variable(3).value;
                              'values4'       {'real','cell' } []     defdes.variable(4).value;
                              'vartype1'      'string' {'categorical' 'continuous'}  defdes.variable(1).vartype;
                              'vartype2'      'string' {'categorical' 'continuous'}  defdes.variable(2).vartype;
                              'vartype3'      'string' {'categorical' 'continuous'}  defdes.variable(3).vartype;
                              'vartype4'      'string' {'categorical' 'continuous'}  defdes.variable(4).vartype;
                              'name'          'string'    {}     defdes.name;
                              'filepath'      'string'    {}     defdes.filepath;
                              'datselect'     'cell'      {}     defdes.include;
                              'dataselect'    'cell'      {}     {};
                              'subjselect'    'cell'      {}     defdes.cases.value;
                              'delfiles'      'string'    { 'on','off', 'limited' } 'off';
                              'verbose'       'string'    { 'on','off' } 'on';
                              'defaultdesign' 'string'    { 'on','off','forceoff'} fastif(nargin < 3, 'on', 'off') }, ...
                              'std_makedesign', 'ignore');
if ischar(opt), error(opt); end
if ~isempty(opt.dataselect), opt.datselect = opt.dataselect; end
if strcmpi(opt.variable1, 'none'), opt.variable1 = ''; end
if strcmpi(opt.variable2, 'none'), opt.variable2 = ''; end
%if iscell(opt.values1), for i = 1:length(opt.values1), if iscell(opt.values1{i}), opt.values1{i} = cell2str(opt.values1{i}); end; end; end
%if iscell(opt.values2), for i = 1:length(opt.values2), if iscell(opt.values2{i}), opt.values2{i} = cell2str(opt.values2{i}); end; end; end
   
% check inputs
% ------------
[indvars, indvarvals, ~, paired] = std_getindvar(STUDY);
if isfield(STUDY.datasetinfo, 'trialinfo')
     alltrialinfo = { STUDY.datasetinfo.trialinfo };
     dattrialselect = cellfun(@(x)([1:length(x)]), alltrialinfo, 'uniformoutput', false);
else alltrialinfo = cell(length(STUDY.datasetinfo));
     for i=1:length(ALLEEG), dattrialselect{i} = [1:ALLEEG(i).trials]; end
end

% get values for each independent variable
% ----------------------------------------
m1 = strmatch(opt.variable1, indvars, 'exact'); if isempty(m1), opt.variable1 = ''; end
m2 = strmatch(opt.variable2, indvars, 'exact'); if isempty(m2), opt.variable2 = ''; end
m3 = strmatch(opt.variable3, indvars, 'exact'); if isempty(m3), opt.variable3 = ''; end
m4 = strmatch(opt.variable4, indvars, 'exact'); if isempty(m4), opt.variable4 = ''; end
if isempty(opt.values1) && ~isempty(opt.variable1), opt.values1 = indvarvals{m1}; end
if isempty(opt.values2) && ~isempty(opt.variable2), opt.values2 = indvarvals{m2}; end
if isempty(opt.values3) && ~isempty(opt.variable3), opt.values3 = indvarvals{m3}; end
if isempty(opt.values4) && ~isempty(opt.variable4), opt.values4 = indvarvals{m4}; end
if isempty(opt.variable1), opt.values1 = { '' }; end
if isempty(opt.variable2), opt.values2 = { '' }; end
if isempty(opt.variable3), opt.values3 = { '' }; end
if isempty(opt.variable4), opt.values4 = { '' }; end

% build command list for history
% ------------------------------
listcom = { 'name' opt.name 'delfiles' opt.delfiles 'defaultdesign' opt.defaultdesign };
if ~isempty(opt.variable1), listcom = { listcom{:} 'variable1' opt.variable1 'values1' opt.values1 'vartype1' opt.vartype1 }; end
if ~isempty(opt.variable2), listcom = { listcom{:} 'variable2' opt.variable2 'values2' opt.values2 'vartype2' opt.vartype2 }; end
if ~isempty(opt.variable3), listcom = { listcom{:} 'variable3' opt.variable3 'values3' opt.values3 'vartype3' opt.vartype3 }; end
if ~isempty(opt.variable4), listcom = { listcom{:} 'variable4' opt.variable4 'values4' opt.values4 'vartype4' opt.vartype4 }; end
if ~isempty(opt.subjselect),  listcom = { listcom{:} 'subjselect'  opt.subjselect }; end
if ~isempty(opt.datselect),   listcom = { listcom{:} 'datselect'  opt.datselect }; end
if ~isempty(opt.filepath),    listcom = { listcom{:} 'filepath'  opt.filepath }; end
if ~isempty(opt.datselect),   listcom = { listcom{:} 'datselect'  opt.datselect }; end
    
% select specific subjects
% ------------------------
datsubjects = { STUDY.datasetinfo.subject };
if ~isempty(opt.subjselect)
    allsubjects = opt.subjselect;
else
    allsubjects = unique_bc( datsubjects );
end
    
% delete design files
% ---------------------
if strcmpi(opt.delfiles, 'on')
    myfprintf(opt.verbose, 'Deleting all files pertaining to design %d\n', designind);
    for index = 1:length(ALLEEG)
        files = fullfile(ALLEEG(index).filepath, sprintf(opt.verbose, 'design%d*.*', designind));
        files = dir(files);
        for indf = 1:length(files)
            delete(fullfile(ALLEEG(index).filepath, files(indf).name));
        end
    end
elseif strcmpi(opt.delfiles, 'limited')
    myfprintf(opt.verbose, 'Deleting all files for STUDY design %d\n', designind);
    for index = 1:length(STUDY.design(designind).cell)
        filedir = [ STUDY.design(designind).cell(index).filebase '.dat*' ];
        filepath = fileparts(filedir);
        files = dir(filedir);
        for indf = 1:length(files)
            %disp(fullfile(filepath, files(indf).name));
            delete(fullfile(filepath, files(indf).name));
        end
    end
    for index = 1:length(STUDY.design(designind).cell)
        filedir = [ STUDY.design(designind).cell(index).filebase '.ica*' ];
        filepath = fileparts(filedir);
        files = dir(filedir);
        for indf = 1:length(files)
            %disp(fullfile(filepath, files(indf).name));
            delete(fullfile(filepath, files(indf).name));
        end
    end
end

% preselect data
% --------------
datselect = [1:length(STUDY.datasetinfo)];
if ~isempty(opt.datselect)
    myfprintf(opt.verbose, 'Data preselection for STUDY design\n');
    for ind = 1:2:length(opt.datselect)
        [ dattmp, dattrialstmp ] = std_selectdataset( STUDY, ALLEEG, opt.datselect{ind}, opt.datselect{ind+1});
        datselect      = intersect_bc(datselect, dattmp);
        dattrialselect = intersectcell(dattrialselect, dattrialstmp);
    end
end
datselect = intersect_bc(datselect, strmatchmult(allsubjects, datsubjects));

des.name              = opt.name;
des.filepath          = opt.filepath;
if ~isempty(opt.variable1)
    des.variable(1).label   = opt.variable1;
    des.variable(1).value   = opt.values1;
    des.variable(1).vartype = opt.vartype1;
    des.variable(1).pairing  = paired{m1};
end
if ~isempty(opt.variable2)
    des.variable(2).label   = opt.variable2;
    des.variable(2).value   = opt.values2;
    des.variable(2).vartype = opt.vartype2;
    des.variable(2).pairing  = paired{m2};
end
if ~isempty(opt.variable3)
    des.variable(3).label   = opt.variable3;
    des.variable(3).value   = opt.values3;
    des.variable(3).vartype = opt.vartype3;
    des.variable(3).pairing  = paired{m3};
end
if ~isempty(opt.variable4)
    des.variable(4).label   = opt.variable4;
    des.variable(4).value   = opt.values4;
    des.variable(4).vartype = opt.vartype4;
    des.variable(4).pairing  = paired{m4};
end
des.include             = opt.datselect;
des.cases.label = 'subject';
des.cases.value = opt.subjselect;
if isempty(des.cases.value) des.cases.value = STUDY.subject; end

fieldorder = { 'name' 'filepath' 'variable' 'cases' 'include' };
if ~isfield(des, 'variable'), des.variable = []; end
des        = orderfields(des, fieldorder);
try
    STUDY.design = orderfields(STUDY.design, fieldorder);
catch
    STUDY.design = [];
end

if ~isfield(STUDY, 'design') || isempty(STUDY.design)
    STUDY.design = des;
else
    STUDY.design(designind) = des;  % fill STUDY.design
end

STUDY = std_addvarlevel(STUDY, designind);

% remove empty variables
if designind <= length(STUDY.design)
    for iVar = length(STUDY.design(designind).variable):-1:1
        if isempty(STUDY.design(designind).variable(iVar).label)
            STUDY.design(designind).variable(iVar) = [];
        end
    end
end

% reorder variable putting single value var at the end
% this allows selecting some specific conditions and still processing the
% data in EEGLAB standard stat functions
if length(STUDY.design(designind).variable) > 2
    if all(cellfun(@(x)isequal(x, 'categorical'), { STUDY.design(designind).variable(3).vartype })) % only categorical
        varLen = cellfun(@(x)length(x), {STUDY.design(designind).variable.value});
        indLen1 = find(varLen == 1);
        STUDY.design(designind).variable = [ STUDY.design(designind).variable(setdiff(1:length(varLen), indLen1)) STUDY.design(designind).variable(indLen1) ];
    end
end

% select the new design in the STUDY output
% -----------------------------------------
STUDY.currentdesign     = designind;
STUDY = std_selectdesign(STUDY, ALLEEG, designind); 
STUDY.cache = [];

% build output command
% --------------------
com = sprintf('STUDY = std_makedesign(STUDY, ALLEEG, %d, %s);', designind, vararg2str( listcom ) );

% ---------------------------------------------------

% return intersection of cell arrays
% ----------------------------------
function res = intersectcell(a, b, c)
if nargin > 2
    res = intersectcell(a, intersectcell(b, c));
else
    for index = 1:length(a)
        res{index} = intersect_bc(a{index}, b{index});
    end
end

% perform multi strmatch
% ----------------------
function res = strmatchmult(a, b)
    res = [];
    for index = 1:length(a)
        res = [ res strmatch(a{index}, b, 'exact')' ];
    end

% remove blanks
% -------------
function res = rmblk(a)
    if iscell(a), a = cell2str(a); end
    if ~ischar(a), a = num2str(a); end
    res = a;
    res(find(res == ' ')) = '_';
    res(find(res == '\')) = '_';    
    res(find(res == '/')) = '_';    
    res(find(res == ':')) = '_';    
    res(find(res == ';')) = '_';    
    res(find(res == '"')) = '_';    
    
function strval = cell2str(vals);
    strval = vals{1};
    for ind = 2:length(vals)
        strval = [ strval ' - ' vals{ind} ];
    end

function tmpfile = checkfilelength(tmpfile)
    if length(tmpfile) > 120
        tmpfile = tmpfile(1:120);
    end

function myfprintf(verbose, varargin)
    if strcmpi(verbose, 'on'), fprintf(varargin{:}); end
