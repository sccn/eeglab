% std_plodtmat() - plot design matrix and info associated with the study for
% each subject. Designed to be used directly (as a callback function of button plot) from pop_studydesign
%
% Usage:
%   >>  
%
% Inputs:
%      usrdat    - Structure who contain the name of the trial's properties
%         factors: Cell array with the name of the trial's properties {1x15 cell}
%      factorvals: Cell array with the values of the properties in 'factors' for each trial
%        factsubj: 
%        subjects: Cell array with the name of the subjects
%     datasetinfo: datasetinfo
%          design: Structure withe the fields {name,filepath,variable,
%          cases, include,cell, deletepreviousfiles} for each variable in
%          the design
%        filepath: Study filepath
%       numerical: 
%       
%      numdesign - Design selected
%    
% Outputs:
%
% See also: 
%   pop_studydesign,
%
% Author: Ramon Martinez-Cancino, SCCN, 2014
%
% Copyright (C) 2014  Ramon Martinez-Cancino, 
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

function std_plotmat(usrdat,numdesign)

% Set var stuff
design = usrdat.design(numdesign);
listsubj = unique(design.cases.value,'sorted'); % assuming cases.value will be alwas the subjects
setappdata(0,'usrdat',usrdat);
setappdata(0,'numdesign',numdesign);

% Creating GUI
% Positions and settings
%--------------------------------------------------------------------------
mainfig_pos        = [.509  .465  .306  .519];
Text1_pos          = [.124  .946  .191  .024];
Text2_pos          = [.080  .208  .432  .029];
Text3_pos          = [.397  .861  .24   .024];
checkbox_sort_pos  = [.422  .901  .242  .027];
axes_pos           = [.118  .306  .817  .513];
listbox1_pos       = [.12   .018  .817  .178];
popup_subject_pos  = [.124  .896  .224  .035];

GUI_FONTSIZE = 9;
COLOR = [.66 .76 1];
figunits = 'Normalized';

% Main fig
%--------------------------------------------------------------------------
handles.mainfig = figure('MenuBar'         ,'none',...
                             'Name'        ,'Design Matrix',...
                             'NumberTitle' ,'off',...
                             'Units'       ,figunits,...
                             'Color'       ,COLOR,...
                             'Position'    ,mainfig_pos);                 
% Text
%--------------------------------------------------------------------------
handles.Text1 = uicontrol('Style','Text');
set(handles.Text1,'String'          ,'Select subject',...
                  'FontSize'        ,GUI_FONTSIZE,...
                  'Units'           ,figunits,...
                  'BackgroundColor' ,COLOR,...
                  'Position'        ,Text1_pos);

handles.Text2 = uicontrol('Style','Text');
set(handles.Text2,'String'          ,'Subject Variables Design',...
                  'FontSize'        ,GUI_FONTSIZE,...
                  'FontWeight'      ,'bold',...
                  'Units'           ,figunits,...
                  'BackgroundColor' ,COLOR,...
                  'Position'        ,Text2_pos); 
              
handles.Text3 = uicontrol('Style','Text');
set(handles.Text3,'String'          ,'Design Matrix',...
                  'FontSize'        ,GUI_FONTSIZE,...
                  'FontWeight'      ,'bold',...
                  'Units'           ,figunits,...
                  'BackgroundColor' ,COLOR,...
                  'Position'        ,Text3_pos);             
% Edit
%--------------------------------------------------------------------------                        
handles.disp_prop = uicontrol('Style','Edit');
set(handles.disp_prop,'String'          ,'Loading..',...
                      'FontSize'        ,GUI_FONTSIZE,...
                      'Units'           ,figunits,...
                      'BackgroundColor' ,COLOR,...
                      'Enable'          ,'inactive',...
                      'Min'             ,1,...
                      'Max'             ,30,...
                      'Position'        ,listbox1_pos);                                                                        
% Axes
%--------------------------------------------------------------------------
handles.axes1 =  axes('unit', 'normalized', 'position', axes_pos);
handles.axes2 =  axes('unit'          ,'normalized', ...
                      'position'      ,axes_pos, ...  
                      'Position'      ,handles.axes1.Position,...
                      'XAxisLocation' ,'top',...
                      'YAxisLocation' ,'right',...
                      'Color'         ,'none',...
                      'YTick'         ,'',...
                      'YTickLabel'    ,'');
% Checkbox
%--------------------------------------------------------------------------
 handles.checkbox_sort = uicontrol('Style'   ,'checkbox');
 set(handles.checkbox_sort,'String'             ,'Sort trials',...
                           'FontSize'        ,GUI_FONTSIZE,...
                           'Units'           ,figunits,...
                           'BackgroundColor' ,COLOR,...
                           'Position'        ,checkbox_sort_pos,...
                           'Value'           ,1',...
                           'Callback'        ,{@callback_popup_subject,handles}); 

callback_popup_subject('', '', handles);

% Popupmenu
%--------------------------------------------------------------------------
handles.popup_subject = uicontrol('Style','Popupmenu');
set(handles.popup_subject,'String'   ,listsubj,...
                          'FontSize' ,GUI_FONTSIZE,...
                          'Units'    ,figunits,...
                          'Position' ,popup_subject_pos,...
                          'callback' ,{@callback_popup_subject,handles});
                      
%--------------------------------------------------------------------------
% Callbacks
%--------------------------------------------------------------------------
function callback_popup_subject(src,eventdata,handles)

numdesign = getappdata(0,'numdesign');
usrdat    = getappdata(0,'usrdat');

design = usrdat.design(numdesign);
listsubj = unique(design.cases.value,'sorted'); % assuming cases.value will be alwas the subjects
if isfield(handles,'popup_subject')
    indtmp = handles.popup_subject.Value;
else
    indtmp = 1;
end

subj = listsubj{indtmp};

%  Retreiving all trials and values for this subject
cell_subjindx     = find(strcmp({design.cell.case},subj));                 % in design.cell
dsetinfo_subjindx = sort(find(strcmp({usrdat.datasetinfo.subject},subj))); % in usrdat.datasetinfo

ntrials = 0;
for i = 1 : length(dsetinfo_subjindx)
    startendindx(i,1) = ntrials + 1;
    ntrials = ntrials + length(usrdat.datasetinfo(i).trialinfo);
    startendindx(i,2) = ntrials;
end

tmpdmat = NaN(ntrials,length(design.variable));
newntrials = 0;
% Filling temporal design matrix tmpdat(loop per factors)
for i = 1 : length(cell_subjindx)
    tmpdset     = design.cell(cell_subjindx(i)).dataset;
    tmptrials   = design.cell(cell_subjindx(i)).trials{:};
    newntrials     = newntrials + (length(tmptrials));
    
    tmpdsetindx = find(tmpdset==dsetinfo_subjindx);
    
    for j = 1 : length(design.variable)
        facval      = cell2mat(design.cell(cell_subjindx(i)).value(j));
        if isnumeric(facval)
            facval_indx = find(facval == cell2mat(design.variable(j).value));
        else
            facval_indx = find(strcmp(facval,design.variable(j).value));
        end
        tmpdmat(tmptrials + startendindx(tmpdsetindx,1)-1,j) = facval_indx;
    end
end

% Removing NANs (Thanks Cyril!!!)
check = find(sum(isnan(tmpdmat),2));
tmpdmat(check,:) = [];

% Checking checbox to sort/unsort
 if isfield(handles, 'checkbox_sort') && get(handles.checkbox_sort, 'Value')
     [tmpdmat,tmp] = sortrows(tmpdmat,[1:size(tmpdmat,2)]);
 end

% Checking checbox to sort/unsort
if isfield(handles, 'checkbox_sort') & handles.checkbox_sort.Value
    [tmpdmat,tmp] = sortrows(tmpdmat,[1:size(tmpdmat,2)]);
end

% Updating the design info
text2display(1) = {['Design Name: ' design.name]}; %Name of the design

% Name and values of Regressors
for i = 1:length(design.variable) % Loop per condition
    allvarstmp = '[';
    for j = 1 : length(design.variable(i).value)
        allvarstmp = [allvarstmp ' ' num2str(design.variable(i).value{j})];
    end
    allvarstmp = [allvarstmp ' ]'];
    
    ncond(i)           = length(design.variable(i).value);
    text2display(2*i)     = {['Regressor ' num2str(i) ' : ' design.variable(i).label]};
    text2display(2*i + 1) = {['Regressor ' num2str(i) ' values :' allvarstmp]};
end

text2display(2*(i+1))     = {['Regressor ' num2str(i+1) ' : Baseline' ]};
text2display(2*(i+1) + 1) = {['Regressor ' num2str(i+1) ' values : [1]' ]};

%  Display values
set(handles.disp_prop,'String',text2display,'HorizontalAlignment', 'left');

% display(['Number of trials missing in design:' num2str(ntrials-newntrials)]);
% fig1 = handles.axes1;
axes(handles.axes1)
imagesc(tmpdmat); colormap(flipud(colormap('gray')));
xlabel('Regressors', 'FontWeight', 'normal', 'FontSize', 9);
ylabel('Trials','FontWeight', 'Normal', 'FontSize', 9);
set(handles.axes1,'XTick',1:size(tmpdmat,2))
set(handles.axes2','XTick', handles.axes1.XTick,...
                   'XTickLabel', {design.variable.label 'Baseline'},...
                   'XLim', handles.axes1.XLim);
drawnow;
