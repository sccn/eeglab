% std_plodtmat() - plot design matrix and info associated with the study for
%                  each subject. Designed to be used directly (as a callback 
%                  function of button plot) from pop_studydesign
% Usage:
%   >>  
%
% Inputs:
%      usrdat      - Structure who contain the name of the trial's properties
%      factors     - Cell array with the name of the trial's properties {1x15 cell}
%      factorvals  -Cell array with the values of the properties in 'factors' for each trial
%      factsubj 
%      subjects    -Cell array with the name of the subjects
%      datasetinfo -datasetinfo
%      design      -Structure withe the fields {name,filepath,variable,
%                   cases, include, deletepreviousfiles} for each variable
%                   in the design
%      filepath    -Study filepath
%      numerical  
%      numdesign   -Design selected
%    
% Outputs:
%
% See also: 
%   pop_studydesign
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

flag.textdisp = 1;
flag.subj     = 1;
setappdata(0,'flag',flag);
handles.figmode = 0; % Change to 1 if info in GUI. 0 -> info in commandline
sortlist = [{' '} {design.variable.label}]';

% Creating GUI
% Positions and settings
%--------------------------------------------------------------------------
mainfig_pos            = [.509  .465  .306  .519];
Text1_pos              = [.120  .946  .191  .024];
Text2_pos              = [.080  .208  .432  .029];
Text3_pos              = [.397  .855  .24   .024];
Text4_pos              = [.680  .946  .201  .024];
Text5_pos              = [.380  .946  .201  .024];
checkbox_sort_pos      = [.422  .896  .242  .035];

popup_subject_pos      = [.124  .896  .224  .035];
disp_click_pos         = [.7    .896  .183  .035];

if handles.figmode
    axes_pos           = [.118  .306  .817  .513];
    listbox1_pos       = [.12   .018  .817  .178];
else
    axes_pos           = [.118  .07   .817  .750];
end

GUI_FONTSIZE = 10;
AXES_FONTSIZE = 9;
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
% Data in mainfig
setappdata(handles.mainfig,'usrdat',usrdat);
setappdata(handles.mainfig,'numdesign',numdesign);
setappdata(handles.mainfig,'flag',flag);

% Text
%--------------------------------------------------------------------------
handles.Text1 = uicontrol('Style','Text');
set(handles.Text1,'String'          ,'Select subject',...
                  'FontSize'        ,GUI_FONTSIZE,...
                  'Units'           ,figunits,...
                  'BackgroundColor' ,COLOR,...
                  'Position'        ,Text1_pos);
if handles.figmode
handles.Text2 = uicontrol('Style','Text');
set(handles.Text2,'String'          ,'Subject Variables Design',...
                  'FontSize'        ,GUI_FONTSIZE,...
                  'FontWeight'      ,'bold',...
                  'Units'           ,figunits,...
                  'BackgroundColor' ,COLOR,...
                  'Position'        ,Text2_pos);
end           
handles.Text3 = uicontrol('Style','Text');
set(handles.Text3,'String'          ,'Design Matrix',...
                  'FontSize'        ,GUI_FONTSIZE,...
                  'FontWeight'      ,'bold',...
                  'Units'           ,figunits,...
                  'BackgroundColor' ,COLOR,...
                  'Position'        ,Text3_pos);

handles.Text4 = uicontrol('Style','Text');
set(handles.Text4,'String'          ,'Value on Click',...
                  'FontSize'        ,GUI_FONTSIZE,...
                  'Units'           ,figunits,...
                  'BackgroundColor' ,COLOR,...
                  'Position'        ,Text4_pos);
              
handles.Text5 = uicontrol('Style','Text');
set(handles.Text5,'String'          ,'Sort by: ',...
                  'FontSize'        ,GUI_FONTSIZE,...
                  'Units'           ,figunits,...
                  'BackgroundColor' ,COLOR,...
                  'Position'        ,Text5_pos);              
              
% Edit
%--------------------------------------------------------------------------  
if handles.figmode
    handles.disp_prop = uicontrol('Style','Edit');
    set(handles.disp_prop,'String'          ,'Loading..',...
                          'FontSize'        ,GUI_FONTSIZE,...
                          'Units'           ,figunits,...
                          'BackgroundColor' ,COLOR,...
                          'Enable'          ,'inactive',...
                          'Min'             ,1,...
                          'Max'             ,30,...
                          'Position'        ,listbox1_pos); 
 end
                  
handles.edit_dispclick = uicontrol('Style','Edit');
set(handles.edit_dispclick,'String'           ,' ',...
                           'FontSize'         ,GUI_FONTSIZE,...
                           'Units'            ,'Normalized',...
                           'enable'           ,'off',...
                           'ForegroundColor'  , [0 0 0],...
                           'BackgroundColor'  ,[1 1 1],...
                           'Position'         ,disp_click_pos);
                         
% Axes
%--------------------------------------------------------------------------
handles.axes1 =  axes('unit', 'normalized', 'position', axes_pos);
handles.axes2 =  axes('unit'          ,'normalized', ...
                      'Position'      ,axes_pos, ...  
                      'XAxisLocation' ,'top',...
                      'YAxisLocation' ,'right',...
                      'Color'         ,'none',...
                      'YTick'         ,[],...
                      'YTickLabel'    ,'');
                      
% Popupmenu
%--------------------------------------------------------------------------   
handles.popup_sort    = uicontrol('Style','Popupmenu');
handles.popup_subject = uicontrol('Style','Popupmenu');

set(handles.popup_sort,'String'   ,sortlist,...
                          'FontSize' ,GUI_FONTSIZE,...
                          'Units'    ,figunits,...
                          'Position' ,checkbox_sort_pos,...
                          'callback' ,{@callback_popup_subject,handles});      
                      
set(handles.popup_subject,'String'   ,listsubj,...
                          'FontSize' ,GUI_FONTSIZE,...
                          'Units'    ,figunits,...
                          'Position' ,popup_subject_pos,...
                          'callback' ,{@callback_popup_subject,handles});
                      
callback_popup_subject('', '', handles);
              
%--------------------------------------------------------------------------
% Callbacks
%--------------------------------------------------------------------------
function callback_popup_subject(src,eventdata,handles)

AXES_FONTSIZE = 11;
numdesign = getappdata(handles.mainfig,'numdesign');
usrdat    = getappdata(handles.mainfig,'usrdat');
set(handles.edit_dispclick, 'String', ' ', ...
                             'ForegroundColor'  , [0 0 0]);
design = usrdat.design(numdesign);
listsubj = unique(design.cases.value,'sorted'); % assuming cases.value will be alwas the subjects
if isfield(handles,'popup_subject')
    indtmp = get(handles.popup_subject, 'Value');
else
    indtmp = 1;
end
flag = getappdata(handles.mainfig,'flag');
inner_subj_flag = 1;
if flag.subj ~= indtmp
    flag.textdisp = 1;
    flag.subj     = indtmp;
    inner_subj_flag = 0; 
end
subj = listsubj{indtmp};

%  Retreiving all trials and values for this subject
dsetinfo_subjindx = sort(find(strcmp({usrdat.datasetinfo.subject},subj)));  % in usrdat.datasetinfo
trialinfo  = std_combtrialinfo(usrdat.datasetinfo, dsetinfo_subjindx);      % Combining trialinfo

ntrials = 0;
for i = 1 : length(dsetinfo_subjindx)
    startendindx(i,1) = ntrials + 1;
    ntrials = ntrials + length(usrdat.datasetinfo(dsetinfo_subjindx(i)).trialinfo);
    startendindx(i,2) = ntrials;
end

% call function to build design matrix
[tmpdmat allLabels] = std_builddesignmat(design, trialinfo, true);

% Removing NANs (Thanks Cyril!!!)
check = find(sum(isnan(tmpdmat),2));
tmpdmat(check,:) = [];

% Checking checbox to sort/unsort
dmatsortindex = 1:size(tmpdmat, 1);
if get(handles.popup_sort, 'Value') ~= 1
    [tmp,dmatsortindex] = sortrows(tmpdmat,get(handles.popup_sort, 'Value')-1);
     if inner_subj_flag, flag.textdisp = 0; end
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
if handles.figmode
    set(handles.disp_prop,'String',text2display,'HorizontalAlignment', 'left','FontSize', AXES_FONTSIZE);
elseif  flag.textdisp 
    display('----------------------------------------')
    display(['--- Subject ' subj ' Variables Design---'])
    disp(text2display(:));
    display('----------------------------------------')
    handles.notdisp = 0;
end
axes(handles.axes1)
handles.figure = imagesc(tmpdmat(dmatsortindex,:)); colormap(flipud(colormap('gray'))); % Normalizing before show (by cols)
xlabel('Regressors',...
       'FontWeight', 'normal',...
       'FontSize', AXES_FONTSIZE);
ylabel('Trials',...
       'FontWeight', 'Normal',...
       'FontSize', AXES_FONTSIZE);
set(handles.axes1,'XTick',1:size(tmpdmat,2))
set(handles.axes2','XTick', get(handles.axes1, 'XTick'),...
                   'XTickLabel', allLabels,...
                   'XLim', get(handles.axes1, 'XLim'),...
                   'FontSize', AXES_FONTSIZE);
set(handles.figure, 'ButtonDownFcn', {@callback_dmatclick,handles,tmpdmat(dmatsortindex,:),design})
setappdata(handles.mainfig,'flag',flag)
drawnow;

function callback_dmatclick(src,eventdata,handles,tmpdmat,design)

axesHandle  = get(src,'Parent');
coordinates = get(axesHandle,'CurrentPoint');
coordinates = round(coordinates(1,1:2));
tmpdmatval = tmpdmat(coordinates(2), coordinates(1));
if  coordinates(1) <= length(design.variable)
    if isfield(design.variable(coordinates(1)),'vartype') && strcmp(design.variable(coordinates(1)).vartype, 'categorical')
        val = design.variable(coordinates(1)).value(tmpdmatval);
    else
        % Note: When an old STUDY is used withouth modifying the variables,
        % the structure 'STUDY.design' does not have the field 'vartype'.
        val = tmpdmatval;
    end
else
    val = 1;
end
if isnumeric(val)
    val = num2str(val);
end
set(handles.edit_dispclick, 'String', val, ...
                             'ForegroundColor'  , [0 0 0]);