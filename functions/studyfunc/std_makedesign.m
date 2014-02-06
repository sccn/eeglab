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
%  'pairing1' - ['on'|'off'] the nature of the 'variable1' contrast. 
%               For example, to compare two conditions recorded from the 
%               same group of 10 subjects, the 'variable1','condition' design 
%               elements are paired ('on') since each dataset for one
%               condition has a corresponding dataset from the same subject 
%               in the second condition. If the two conditions were recorded 
%               from different groups of subjects, the variable1 'condition' 
%               would be unpaired ('off') {default: 'on'}
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
%  'pairing2' - ['on'|'off'] type of statistics for variable2 
%               (default: 'on'}
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

function [STUDY com] = std_makedesign(STUDY, ALLEEG, designind, varargin)

if nargin < 2
    help std_makedesign;
    return;
end;
if nargin < 3 
    designind = 1;
end;

defdes.name = sprintf('STUDY.design %d', designind);
defdes.cases.label = 'subject';
defdes.cases.value = {};
defdes.variable(1).label = 'condition';
defdes.variable(2).label = 'group';
defdes.variable(1).value = {};
defdes.variable(2).value = {};
defdes.variable(1).pairing = 'on';
defdes.variable(2).pairing = 'on';
defdes.filepath = '';
defdes.include = {};
orivarargin = varargin;
if ~isempty(varargin) && isstruct(varargin{1})
    defdes = varargin{1};
    varargin(1) = [];
end;
if isempty(defdes.filepath), defdes.filepath = ''; end;
opt = finputcheck(varargin,  {'variable1'     'string'    []     defdes.variable(1).label;
                              'variable2'     'string'    []     defdes.variable(2).label;
                              'values1'       {'real','cell' } []     defdes.variable(1).value;
                              'values2'       {'real','cell' } []     defdes.variable(2).value;
                              'pairing1'      'string'    []     defdes.variable(1).pairing;
                              'pairing2'      'string'    []     defdes.variable(2).pairing;
                              'name'          'string'    {}     defdes.name;
                              'filepath'      'string'    {}     defdes.filepath;
                              'datselect'     'cell'      {}     defdes.include;
                              'dataselect'    'cell'      {}     {};
                              'subjselect'    'cell'      {}     defdes.cases.value;
                              'delfiles'      'string'    { 'on','off', 'limited' } 'off';
                              'verbose'       'string'    { 'on','off' } 'on';
                              'defaultdesign' 'string'    { 'on','off','forceoff'} fastif(nargin < 3, 'on', 'off') }, ...
                              'std_makedesign');
if isstr(opt), error(opt); end;
if ~isempty(opt.dataselect), opt.datselect = opt.dataselect; end;
if strcmpi(opt.variable1, 'none'), opt.variable1 = ''; end;
if strcmpi(opt.variable2, 'none'), opt.variable2 = ''; end;
%if iscell(opt.values1), for i = 1:length(opt.values1), if iscell(opt.values1{i}), opt.values1{i} = cell2str(opt.values1{i}); end; end; end;
%if iscell(opt.values2), for i = 1:length(opt.values2), if iscell(opt.values2{i}), opt.values2{i} = cell2str(opt.values2{i}); end; end; end;
    
% build command list for history
% ------------------------------
listcom = { 'variable1' opt.variable1 'variable2' opt.variable2 'name' opt.name 'pairing1' opt.pairing1 'pairing2' opt.pairing2 'delfiles' opt.delfiles 'defaultdesign' opt.defaultdesign };
if ~isempty(opt.values1), listcom = { listcom{:} 'values1' opt.values1 }; end;
if ~isempty(opt.values2), listcom = { listcom{:} 'values2' opt.values2 }; end;
if ~isempty(opt.subjselect),  listcom = { listcom{:} 'subjselect'  opt.subjselect }; end;
if ~isempty(opt.datselect),   listcom = { listcom{:} 'datselect'  opt.datselect }; end;
if ~isempty(opt.filepath),    listcom = { listcom{:} 'filepath'  opt.filepath }; end;
if ~isempty(opt.datselect),   listcom = { listcom{:} 'datselect'  opt.datselect }; end;
    
% select specific subjects
% ------------------------
datsubjects = { STUDY.datasetinfo.subject };
if ~isempty(opt.subjselect)
    allsubjects = opt.subjselect;
else
    allsubjects = unique_bc( datsubjects );
end;
    
% delete design files
% ---------------------
if strcmpi(opt.delfiles, 'on')
    myfprintf(opt.verbose, 'Deleting all files pertaining to design %d\n', designind);
    for index = 1:length(ALLEEG)
        files = fullfile(ALLEEG(index).filepath, sprintf(opt.verbose, 'design%d*.*', designind));
        files = dir(files);
        for indf = 1:length(files)
            delete(fullfile(ALLEEG(index).filepath, files(indf).name));
        end;
    end;
elseif strcmpi(opt.delfiles, 'limited')
    myfprintf(opt.verbose, 'Deleting all files for STUDY design %d\n', designind);
    for index = 1:length(STUDY.design(designind).cell)
        filedir = [ STUDY.design(designind).cell(index).filebase '.dat*' ];
        filepath = fileparts(filedir);
        files = dir(filedir);
        for indf = 1:length(files)
            %disp(fullfile(filepath, files(indf).name));
            delete(fullfile(filepath, files(indf).name));
        end;
    end;
    for index = 1:length(STUDY.design(designind).cell)
        filedir = [ STUDY.design(designind).cell(index).filebase '.ica*' ];
        filepath = fileparts(filedir);
        files = dir(filedir);
        for indf = 1:length(files)
            %disp(fullfile(filepath, files(indf).name));
            delete(fullfile(filepath, files(indf).name));
        end;
    end;
end;

% check inputs
% ------------
[indvars indvarvals ] = std_getindvar(STUDY);
if isfield(STUDY.datasetinfo, 'trialinfo')
     alltrialinfo = { STUDY.datasetinfo.trialinfo };
     dattrialselect = cellfun(@(x)([1:length(x)]), alltrialinfo, 'uniformoutput', false);
else alltrialinfo = cell(length(STUDY.datasetinfo));
     for i=1:length(ALLEEG), dattrialselect{i} = [1:ALLEEG(i).trials]; end;
end;

% get values for each independent variable
% ----------------------------------------
m1 = strmatch(opt.variable1, indvars, 'exact'); if isempty(m1), opt.variable1 = ''; end;
m2 = strmatch(opt.variable2, indvars, 'exact'); if isempty(m2), opt.variable2 = ''; end;
if isempty(opt.values1) && ~isempty(opt.variable1), opt.values1 = indvarvals{m1}; end;
if isempty(opt.values2) && ~isempty(opt.variable2), opt.values2 = indvarvals{m2}; end;
if isempty(opt.variable1), opt.values1 = { '' }; end;
if isempty(opt.variable2), opt.values2 = { '' }; end;

% preselect data
% --------------
datselect = [1:length(STUDY.datasetinfo)];
if ~isempty(opt.datselect)
    myfprintf(opt.verbose, 'Data preselection for STUDY design\n');
    for ind = 1:2:length(opt.datselect)
        [ dattmp dattrialstmp ] = std_selectdataset( STUDY, ALLEEG, opt.datselect{ind}, opt.datselect{ind+1});
        datselect      = intersect_bc(datselect, dattmp);
        dattrialselect = intersectcell(dattrialselect, dattrialstmp);
    end;
end;
datselect = intersect_bc(datselect, strmatchmult(allsubjects, datsubjects));

% get the dataset and trials for each of the ind. variable
% --------------------------------------------------------
ns  = length(allsubjects);
nf1 = max(1,length(opt.values1));
nf2 = max(1,length(opt.values2));
myfprintf(opt.verbose, 'Building STUDY design\n');
for n1 = 1:nf1, [ dats1{n1} dattrials1{n1} ] = std_selectdataset( STUDY, ALLEEG, opt.variable1, opt.values1{n1}, fastif(strcmpi(opt.verbose, 'on'), 'verbose', 'silent')); end;
for n2 = 1:nf2, [ dats2{n2} dattrials2{n2} ] = std_selectdataset( STUDY, ALLEEG, opt.variable2, opt.values2{n2}, fastif(strcmpi(opt.verbose, 'on'), 'verbose', 'silent')); end;

% detect files from old format
% ----------------------------
if ~strcmpi(opt.defaultdesign, 'forceoff') && isempty(opt.filepath)
    if designind == 1
        if strcmpi(opt.defaultdesign, 'off')
            if isfield(STUDY, 'design') && ( ~isfield(STUDY.design, 'cell') || ~isfield(STUDY.design(1).cell, 'filebase') )
                 opt.defaultdesign = 'on';
            end;
        end;
        if isempty(dir(fullfile(ALLEEG(1).filepath, [ ALLEEG(1).filename(1:end-4) '.dat*' ]))) && ...
                isempty(dir(fullfile(ALLEEG(1).filepath, [ ALLEEG(1).filename(1:end-4) '.ica*' ])))
            opt.defaultdesign = 'off';
        end;
    else
        opt.defaultdesign = 'off';
    end;
else
    opt.defaultdesign = 'off';
end;
    
% scan subjects and conditions
% ----------------------------
count = 1;
for n1 = 1:nf1
    for n2 = 1:nf2
        % create design for this set of conditions and subject
        % ----------------------------------------------------
        datasets = intersect_bc(intersect(dats1{n1}, dats2{n2}), datselect);
        if ~isempty(datasets)
            subjects = unique_bc(datsubjects(datasets));
            for s = 1:length(subjects)
                datsubj = datasets(strmatch(subjects{s}, datsubjects(datasets), 'exact'));
                des.cell(count).dataset   = datsubj;
                des.cell(count).trials    = intersectcell(dattrialselect(datsubj), dattrials1{n1}(datsubj), dattrials2{n2}(datsubj));
                des.cell(count).value     = { opt.values1{n1} opt.values2{n2} };
                des.cell(count).case      = subjects{s};
                defaultFile = fullfile(ALLEEG(datsubj(1)).filepath, ALLEEG(datsubj(1)).filename(1:end-4));
                dirres1 = dir( [ defaultFile '.dat*' ] );
                dirres2 = dir( [ defaultFile '.ica*' ] );
                if strcmpi(opt.defaultdesign, 'on') && (~isempty(dirres1) || ~isempty(dirres2)) && isempty(opt.filepath)
                     des.cell(count).filebase = defaultFile;
                else
                    if isempty(rmblk(opt.values1{n1})),    txtval = rmblk(opt.values2{n2});
                    elseif isempty(rmblk(opt.values2{n2})) txtval = rmblk(opt.values1{n1});
                    else txtval =  [ rmblk(opt.values1{n1}) '_' rmblk(opt.values2{n2}) ];
                    end;
                    if ~isempty(txtval),      txtval = [ '_' txtval ]; end;
                    if ~isempty(subjects{s}), txtval = [ '_' rmblk(subjects{s}) txtval ]; end;
                    if isempty(opt.filepath), tmpfilepath = ALLEEG(datsubj(1)).filepath; else tmpfilepath = opt.filepath; end;
                    des.cell(count).filebase = fullfile(tmpfilepath, [ 'design' int2str(designind) txtval ] );
                    des.cell(count).filebase = checkfilelength(des.cell(count).filebase);
                end;
                count = count+1;
            end;
        end;
    end;
end;

% create other fields for the design
% ----------------------------------
if exist('des') ~= 1
    error( [ 'One of your design is empty. This could be because the datasets/subjects/trials' 10 ...
             'you have selected do not contain any of the selected independent variables values.' 10 ...
             'Check your data and datasets carefully for any missing information.' ]);
else    
    % check for duplicate entries in filebase
    % ---------------------------------------
    if length( { des.cell.filebase } ) > length(unique({ des.cell.filebase }))
        if ~isempty(findstr('design_', des.cell(1).filebase))
            error('There is a problem with your STUDY, contact EEGLAB support');
        else
            disp('Duplicate entry detected in new design, reinitializing design with new file names');
            [STUDY com] = std_makedesign(STUDY, ALLEEG, designind, orivarargin{:}, 'defaultdesign', 'forceoff');
            return;
        end
    end;

    %allval1 = unique_bc(cellfun(@(x)x{1}, { des.cell.value }, 'uniformoutput', false));
    %allval2 = unique_bc(cellfun(@(x)x{2}, { des.cell.value }, 'uniformoutput', false));
    des.name              = opt.name;
    des.filepath          = opt.filepath;
    des.variable(1).label = opt.variable1;
    des.variable(2).label = opt.variable2;
    des.variable(1).pairing = opt.pairing1;
    des.variable(2).pairing = opt.pairing2;
    des.variable(1).value   = opt.values1;
    des.variable(2).value   = opt.values2;
    des.include             = opt.datselect;
    des.cases.label = 'subject';
    des.cases.value = unique_bc( { des.cell.case });
end;

fieldorder = { 'name' 'filepath' 'variable' 'cases' 'include' 'cell' };
des        = orderfields(des, fieldorder);
try, 
    STUDY.design = orderfields(STUDY.design, fieldorder);
catch,
    STUDY.design = [];
end;

if ~isfield(STUDY, 'design') || isempty(STUDY.design)
    STUDY.design = des;
else
    STUDY.design(designind) = des;  % fill STUDY.design
end;

% select the new design in the STUDY output
% -----------------------------------------
STUDY.currentdesign     = designind;
STUDY = std_selectdesign(STUDY, ALLEEG, designind); 

% build output command
% --------------------
com = sprintf('STUDY = std_makedesign(STUDY, ALLEEG, %d, %s);', designind, vararg2str( listcom ) );

% ---------------------------------------------------

% return intersection of cell arrays
% ----------------------------------
function res = intersectcell(a, b, c);
if nargin > 2
    res = intersectcell(a, intersectcell(b, c));
else
    for index = 1:length(a)
        res{index} = intersect_bc(a{index}, b{index});
    end;
end;

% perform multi strmatch
% ----------------------
function res = strmatchmult(a, b);
    res = [];
    for index = 1:length(a)
        res = [ res strmatch(a{index}, b, 'exact')' ];
    end;

% remove blanks
% -------------
function res = rmblk(a);
    if iscell(a), a = cell2str(a); end;
    if ~isstr(a), a = num2str(a); end;
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
    end;

function tmpfile = checkfilelength(tmpfile);
    if length(tmpfile) > 120,
        tmpfile = tmpfile(1:120);
    end;

function myfprintf(verbose, varargin);
    if strcmpi(verbose, 'on'), fprintf(varargin{:}); end;