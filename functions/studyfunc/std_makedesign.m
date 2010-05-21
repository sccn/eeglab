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
%   designind   - [integer > 0] index (number) of the new STUDY design {default: 1}
%
% Optional inputs:
%  'name'     - ['string'] mnemonic design name (ex: 'Targets only') 
%               {default: 'Design d', where d = designind} 
%  'indvar1'  - ['string'] - first independent variable or contrast. Must 
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
%  'stat1'    - ['paired'|'unpaired'] the nature of the 'indvar1' contrast. 
%               For example, to compare two conditions recorded from the 
%               same group of 10 subjects, the 'indvar1','condition' design 
%               elements are 'paired' since each dataset for one
%               condition has a corresponding dataset from the same subject 
%               in the second condition. If the two conditions were recorded 
%               from different groups of subjects, the indvar1 'condition' 
%               would be 'unpaired' {default: 'paired'}
%  'indval1'  - {cell array of 'strings'} - 'indvar1' instances to include 
%               in the design. For example, if 'indvar1' is 'condition'and 
%               three values for 'condition' (e.g., 'a' , 'b', and 'c')
%               are listed in STUDY.datasetinfo, then 'indval1', { 'a' 'b' } 
%               will contrast conditions 'a' and 'b', and datasets for 
%               condition 'c' will be ignored. To combine conditions, use 
%               nested '{}'s. For example, to combine conditions 'a' and 
%               'b' into one condition and contrast it to condition 'c', 
%               specify  'indval1', { { 'a' 'b' } 'c' } {default: all values 
%               of 'indvar1' in STUDY.datasetinfo}
%  'indvar2'  - ['string'] - second independent variable name, if any. Typically, 
%               this might refer to ('unpaired') subject group or (typically 
%               'paired') session number, etc.
%  'stat2'    - ['paired'|'unpaired'] type of statistics for indvar2 
%               (default: 'paired'}
%  'indval2'  - {cell array of 'strings'} - indvar2 values to include in the 
%               design {default: all}. Here, 'var[12]' must be field names 
%               in STUDY.datasetinfo or  STUDY.datasetinfo.trialinfo. 
%  'datselect'  - {cell array} select specific datasets and/or trials: 'datselect',
%               {'var1' {'vals'}  'var2' {'vas'}}. Selected datasets must 
%               meet all the specified conditions. For example, 'datselect',  
%               { 'condition' { 'a' 'b' } 'group' { 'g1' 'g2' } }  will 
%               select only datasets from conditions 'a' OR 'b' AND only 
%               subjects in groups 'g1' OR 'g2'. If 'subjselect' is also 
%               specified, only datasets meeting both criteria are included. 
%               'indvar1' and 'indvar2' will only consider
%               the values after they have passed through 'datselect' and 
%               'subjselect'. For instance, if conditions { 'a' 'b' 'c' } 
%               exist and conditions 'a' is removed by 'datselect', the only 
%               two conditions that will be considered are 'b' and 'c' 
%               (which is then equivalent to using 'indvar1vals' to specify
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
%
% Output:
%      STUDY - The input STUDY with a new design added and designated 
%              as the current design.
%
% Examples:
%  STUDY = std_makedesign(STUDY, ALLEEG); % make default design
%
%  % create design with 1 independent variable equal to 'condition'
%  STUDY = std_makedesign(STUDY, ALLEEG, 2, 'indvar1', 'condition');
%
%  % create design with 1 independent variable equal to 'condition'
%  % but only consider the sub-condition 'stim1' and 'stim2' - of course
%  % such conditions must be present in the STUDY
%  STUDY = std_makedesign(STUDY, ALLEEG, 2, 'indvar1', 'condition', ...
%                    'indval1', { 'stim1' 'stim2' }); 
%
%  % create design and only include subject 's1' and 's2'
%  STUDY = std_makedesign(STUDY, ALLEEG, 2, 'indvar1', 'condition', ...
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
defdes.condition   = {};
defdes.group       = {};
defdes.subject     = {};
defdes.indvar1     = 'condition';
defdes.indvar2     = 'group';
defdes.statvar1    = 'paired';
defdes.statvar2    = 'paired';
defdes.includevarlist = {};
if ~isempty(varargin) && isstruct(varargin{1})
    defdes = varargin{1};
    varargin(1) = [];
end;
opt = finputcheck(varargin,  {'indvar1'       'string'    []     defdes.indvar1;
                              'indvar2'       'string'    []     defdes.indvar2;
                              'indval1'       'cell'      []     defdes.condition;
                              'indval2'       'cell'      []     defdes.group;
                              'stat1'         'string'    []     defdes.statvar1;
                              'stat2'         'string'    []     defdes.statvar2;
                              'name'          'string'    {}     defdes.name;
                              'datselect'     'cell'      {}     defdes.includevarlist;
                              'subjselect'    'cell'      {}     defdes.subject;
                              'delfiles'      'string'    { 'on' 'off' } 'off';
                              'defaultdesign' 'string'    { 'on' 'off' } fastif(nargin < 3, 'on', 'off') }, ...
                              'std_makedesign');
if isstr(opt), error(opt); end;
if strcmpi(opt.indvar1, 'none'), opt.indvar1 = ''; end;
if strcmpi(opt.indvar2, 'none'), opt.indvar2 = ''; end;
for i = 1:length(opt.indval1), if iscell(opt.indval1{i}), opt.indval1{i} = cell2str(opt.indval1{i}); end; end;
for i = 1:length(opt.indval2), if iscell(opt.indval2{i}), opt.indval2{i} = cell2str(opt.indval2{i}); end; end;
    
% build command list for history
% ------------------------------
listcom = { 'indvar1' opt.indvar1 'indvar2' opt.indvar2 'name' opt.name };
if ~isempty(opt.indval1), listcom = { listcom{:} 'indvarvals1' opt.indval1 }; end;
if ~isempty(opt.indval2), listcom = { listcom{:} 'indvarvals2' opt.indval2 }; end;
if ~isempty(opt.subjselect),  listcom = { listcom{:} 'subjselect'  opt.subjselect }; end;
if ~isempty(opt.datselect),   listcom = { listcom{:} 'dataselect'  opt.datselect }; end;
    
% select specific subjects
% ------------------------
datsubjects = { STUDY.datasetinfo.subject };
if ~isempty(opt.subjselect)
    allsubjects = opt.subjselect;
else
    allsubjects = unique( datsubjects );
end;
    
% delete design files
% ---------------------
if strcmpi(opt.delfiles, 'on')
    fprintf('Deleting all files for STUDY design %d\n', designind);
    for index = 1:length(ALLEEG)
        files = fullfile(ALLEEG(index).filepath, sprintf('design%d*.*', designind));
        files = dir(files);
        for indf = 1:length(files)
            delete(fullfile(ALLEEG(index).filepath, files(indf).name));
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
m1 = strmatch(opt.indvar1, indvars, 'exact'); if isempty(m1), opt.indvar1 = ''; end;
m2 = strmatch(opt.indvar2, indvars, 'exact'); if isempty(m2), opt.indvar2 = ''; end;
if isempty(opt.indval1) && ~isempty(opt.indvar1), opt.indval1 = indvarvals{m1}; end;
if isempty(opt.indval2) && ~isempty(opt.indvar2), opt.indval2 = indvarvals{m2}; end;
if isempty(opt.indvar1), opt.indval1 = { '' }; end;
if isempty(opt.indvar2), opt.indval2 = { '' }; end;

% preselect data
% --------------
datselect      = [1:length(STUDY.datasetinfo)];
if ~isempty(opt.datselect)
    fprintf('Data preselection for STUDY design\n');
    for ind = 1:2:length(opt.datselect)
        [ dattmp dattrialstmp ] = std_selectdataset( STUDY, opt.datselect{ind}, opt.datselect{ind+1});
        datselect      = intersect(datselect, dattmp);
        dattrialselect = intersectcell(dattrialselect, dattrialstmp);
    end;
end;
datselect = intersect(datselect, strmatchmult(allsubjects, datsubjects));

% get the dataset and trials for each of the ind. variable
% --------------------------------------------------------
ns  = length(allsubjects);
nf1 = max(1,length(opt.indval1));
nf2 = max(1,length(opt.indval2));
fprintf('Building STUDY design\n');
for n1 = 1:nf1, [ dats1{n1} dattrials1{n1} ] = std_selectdataset( STUDY, ALLEEG, opt.indvar1, opt.indval1{n1}); end;
for n2 = 1:nf2, [ dats2{n2} dattrials2{n2} ] = std_selectdataset( STUDY, ALLEEG, opt.indvar2, opt.indval2{n2}); end;

% detect files from old format
% ----------------------------
if strcmpi(opt.defaultdesign, 'on')
    if isempty(dir(fullfile(ALLEEG(1).filepath, '*.dat*'))) || isempty(dir(fullfile(ALLEEG(1).filepath, '*.ica*')))
        opt.defaultdesign = 'off';
    end;
end;

% scan subjects and conditions
% ----------------------------
count = 1;
for n1 = 1:nf1
    for n2 = 1:nf2
        % create design for this set of conditions and subject
        % ----------------------------------------------------
        datasets = intersect(intersect(dats1{n1}, dats2{n2}), datselect);
        if ~isempty(datasets)
            subjects = unique(datsubjects(datasets));
            for s = 1:length(subjects)
                datsubj = datasets(strmatch(subjects{s}, datsubjects(datasets), 'exact'));
                des.setinfo(count).setindex  = datsubj;
                des.setinfo(count).condition = opt.indval1{n1};
                des.setinfo(count).group     = opt.indval2{n2};
                des.setinfo(count).subject   = subjects{s};
                if strcmpi(opt.defaultdesign, 'on')
                     des.setinfo(count).filebase = fullfile(ALLEEG(datsubj(1)).filepath, ALLEEG(datsubj(1)).filename(1:end-4));
                else des.setinfo(count).filebase = fullfile(ALLEEG(datsubj(1)).filepath, ...
                      [ 'design' int2str(designind) '_' subjects{s} '_' rmblk(opt.indval1{n1}) '_' rmblk(opt.indval2{n2}) ] );
                end;
                des.setinfo(count).trialindices = intersectcell(dattrialselect(datsubj), dattrials1{n1}(datsubj), dattrials2{n2}(datsubj));
                count = count+1;
            end;
        end;
    end;
end;

% create other fields for the design
% ----------------------------------
if exist('des') ~= 1
    des = defdes;
    des.name = '';
    disp('Warning: STUDY.design is empty');
else
    des.condition = unique( { des.setinfo.condition });
    des.group     = unique( { des.setinfo.group });
    des.subject   = unique( { des.setinfo.subject });
    if isempty(des.group{1})    , des.group     = {}; end;
    if isempty(des.condition{1}), des.condition = {}; end;
    des.name      = opt.name;
    des.indvar1   = opt.indvar1;
    des.indvar2   = opt.indvar2;
    des.statvar1  = opt.stat1;
    des.statvar2  = opt.stat2;
    des.includevarlist = opt.datselect;
end;

fieldorder = { 'name' 'indvar1' 'indvar2' 'condition' 'group' 'statvar1' 'statvar2' 'subject' 'includevarlist' 'setinfo' };
des          = orderfields(des, fieldorder);
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
if strcmpi(opt.defaultdesign, 'on')
    com = sprintf('STUDY = std_makedesign(STUDY, ALLEEG);' );
else
    com = sprintf('STUDY = std_makedesign(STUDY, ALLEEG, %d, %s);', designind, vararg2str( listcom ) );
end;

% ---------------------------------------------------

% return intersection of cell arrays
% ----------------------------------
function res = intersectcell(a, b, c);
if nargin > 2
    res = intersectcell(a, intersectcell(b, c));
else
    for index = 1:length(a)
        res{index} = intersect(a{index}, b{index});
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
    res = a;
    res(find(res == ' ')) = [];
    
function strval = cell2str(vals);
    strval = vals{1};
    for ind = 2:length(vals)
        strval = [ strval ' - ' vals{ind} ];
    end;
