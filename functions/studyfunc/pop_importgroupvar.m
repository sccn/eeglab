% pop_importgroupvar() - import a group variable from a ASCII file and
% manual entry.
%
% Usage:
%   >>  
%
% Inputs:
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

function [usrdatmod, handles]= pop_importgroupvar(usrdat,numdesign)

if nargin < 2
    help pop_importgroupvar;
    return;
end

% Set var stuff
design = usrdat.design(numdesign);
listsubj = unique(design.cases.value,'sorted'); % assuming cases.value will be alwas the subjects

setappdata(0,'usrdat',usrdat);
setappdata(0,'numdesign',numdesign);
% setappdata(0,'opts',g);

usrdatmod = usrdat;

[filename, pathname, filterindex] = uigetfile({  '*.txt','Text files (*.txt)'; ...
                                                 '*.mat','MAT-files (*.mat)'; ...
                                                 '*.*',  'All Files (*.*)'}, ...
                                                 'Pick a file', ...
                                                 'MultiSelect', 'off');
data = importdata([pathname filename]);
 
% Checking dimensions
if length(data) ~= length(listsubj)
    error('pop_importgroupvar() error: Number of elements in file are inconsistent');
end

% Creating GUI
% Positions and settingss
%--------------------------------------------------------------------------
mainfig_pos           = [.565  .282  .146  .417];  
text_varname_pos      = [.048  .927  .324  .035];
edit_varname_pos      = [.462  .912  .333  .059];
uitable_import_pos    = [.048  .021  .514  .787 ];

button_cancel_pos     = [.59   .024  .357   .056];
button_add_pos        = [.586  .093  .362   .056];

text_vartype_pos      = [.048  .853  .362  .035];
popup_vartype_pos     = [.438  .827  .548  .072];

GUI_FONTSIZE = 9;
COLOR = [.66 .76 1];
figunits = 'Normalized';

% Main fig
%--------------------------------------------------------------------------
handles.mainfig = figure('MenuBar'         ,'none',...
                             'Name'        ,'Import Variable',...
                             'NumberTitle' ,'off',...
                             'Units'       ,figunits,...
                             'Color'       ,COLOR,...
                             'Position'    ,mainfig_pos); 
                         
                                                
% Text
% ......................................................................... 
handles.text_vartype = uicontrol('Style'       ,'text',...
                                  'String'      ,'Variable Type :',...
                                  'FontSize'    ,GUI_FONTSIZE,...
                                  'Units'       ,figunits,...
                                  'BackgroundColor'  ,COLOR,...
                                  'Position'    ,text_vartype_pos);
                             
handles.text_varname= uicontrol('Style'       ,'text',...
                                  'String'      ,'Variable Name : :',...
                                  'FontSize'    ,GUI_FONTSIZE,...
                                  'Units'       ,figunits,...
                                  'BackgroundColor'  ,COLOR,...
                                  'Position'    ,text_varname_pos);   
                              
% Edit
% .........................................................................
handles.edit_varname = uicontrol('Style'            ,'Edit',...
                                 'FontSize'         ,GUI_FONTSIZE,...
                                 'Units'            ,figunits,...
                                 'enable'           ,'on',...
                                 'Position'         ,edit_varname_pos);  
                          
% Uitable
% .........................................................................     
handles.uitable = uitable('units'          ,figunits,...
                          'Position'       ,uitable_import_pos,...
                          'Data'           , data,...
                          'RowName'        ,usrdat.subjects,...
                          'ColumnName'     ,{'Value'},...
                          'ColumnEditable' ,logical(repmat(1,1,length(usrdat.subjects))));
                      
% Popupmenu
% .........................................................................       
handles.popup =uicontrol('Style','popupmenu',...
                         'units'    ,figunits,...
                         'String',{'Categorical';'Continous'},...
                         'Position',popup_vartype_pos);           
                     
% Buttoms
% .........................................................................
handles.button_cancel = uicontrol('Style'       ,'PushButton',...
                                  'String'      ,'Cancel',...
                                  'FontSize'    ,GUI_FONTSIZE,...
                                  'Units'       ,figunits,...
                                  'Position'    ,button_cancel_pos,...
                                  'CallBack'    , {@callback_buton_cancel,handles.mainfig});
                              
handles.button_add     = uicontrol('Style'        ,'PushButton',...
                                   'String'       ,'Add Variable',...
                                   'FontSize'     ,GUI_FONTSIZE,...
                                   'Units'        ,figunits,...
                                   'Position'     ,button_add_pos,...
                                   'CallBack'     ,{@callback_button_add,handles});                     
% .........................................................................

uiwait(handles.mainfig);

% OUTPUTS


%--------------------------------------------------------------------------
% Callbacks
%--------------------------------------------------------------------------  
function callback_buton_cancel(src,eventdata,h)
close(h);

function callback_button_add(src,eventdata,handles)

numdesign = getappdata(0,'numdesign');
usrdat    = getappdata(0,'usrdat');

% retreiving name of variable
varname   = get(handles.edit_varname,'String') ;
if isempty(varname)
    error('pop_importgroupvar() error: Variable name must be provided');
end

% populate fields in STUDY.design
tmpstruct = usrdat.design(numdesign);
if isfield(tmpstruct, 'groupvar')
    if length([tmpstruct.groupvar]) > 1
        
        tmpindx2fill = length([tmpstruct.groupvar]) + 1;
        
    elseif length([tmpstruct.groupvar]) == 1 &...
            all(~isempty(tmpstruct.groupvar.label),...
            ~isempty(tmpstruct.groupvar.value),...
            ~isempty(tmpstruct.groupvar.continuous))
        
        tmpindx2fill = length([tmpstruct.groupvar]) + 1;
        
    elseif length([tmpstruct.groupvar]) == 1 &...
            any(isempty(tmpstruct.groupvar.label),...
            isempty(tmpstruct.groupvar.value),...
            isempty(tmpstruct.groupvar.continuous))
        
        tmpindx2fill = length([tmpstruct.groupvar]);
        
    end
else
    tmpindx2fill = 1;
end

cont = true;    
if get(handles.popup,'Value') ==  1
    cont = false;
end
 
usrdat.design(numdesign).groupvar(tmpindx2fill).label      = varname;                     % design.groupvar.label
usrdat.design(numdesign).groupvar(tmpindx2fill).continuous = cont;                        % design.groupvar.continuous
usrdat.design(numdesign).groupvar(tmpindx2fill).value      = get(handles.uitable,'Data'); % design.groupvar.value

% % populate fields in STUDY.datasetinfo
% usrdat.datasetinfo.(handles.edit_varname.String) = ''
% usrdat.datasetinfo = setfield(usrdat.datasetinfo,handles.edit_varname.String,get(handles.uitable,'Data'));
% 
setappdata(0,'usrdat',usrdat);
close(handles.mainfig);

