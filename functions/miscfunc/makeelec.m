% makeelec() - subroutine to make electrode file in eegplot()
%
% Usage:   >> makeelec(chans) 
%          >> [channames] = makeelec(chans)
%
% Inputs: 
%   chans - number of channels
%
% Author: Colin Humprhies, CNL / Salk Institute, 1996
%
% See also: eegplot()

% Copyright (C) Colin Humphries, CNL / Salk Institute, Aug, 1996
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

% 01-25-02 reformated help & license, added links -ad 

function channames = makeelec(chans);

global is_done_pressed
is_done_pressed = 0;

lines = ceil(chans/2);

fh = 180 + chans*12;
fl = 400;

figure('Units','Pixels','Position',[400 200 fl fh])

fighandle = gcf;

axes('Units','Pixels','Position',[0 0 fl fh],'Visible','off','Xlim',[0 fl],'Ylim',[0 fh]);

text(fl/2,fh-20,'Make Electrode-File','HorizontalAlignment','Center','Fontsize',14)

titleb = uicontrol('Units','Pixels','Style','Edit','Position',...
                     [fl/2 fh-60 .4*fl 20],'UserData',0,'HorizontalAlignment','left');
text(fl/2-20,fh-50,'Filename:','HorizontalAlignment','Right')

text(fl/2,fh-75,'Electrodes','HorizontalAlignment','center')
%line([fl*.2 fl*.8],[fh-85 fh-85],'color','w')

for i = 1:chans

   k = 2*ceil(i/2);

   if (i==k)
      ede(i,:) = uicontrol('Units','Pixels','Style','Edit',...
'Position',[fl*.65 fh-90-k*12 fl*.15 20],'UserData',i,'HorizontalAlignment','left');
      text(fl*.65-15,fh-90-k*12+10,num2str(i),'HorizontalAlignment','right')
   else
      ede(i,:) = uicontrol('Units','Pixels','Style','Edit',...
'Position',[fl*.2 fh-90-k*12 fl*.15 20],'UserData',i,'HorizontalAlignment','left');
      text(fl*.2-15,fh-90-k*12+10,num2str(i),'HorizontalAlignment','right')
   end

end

TIMESTRING = ['chans1973 = get(gco,''UserData'');','OH1973 = findobj(''Style'',''edit'',''UserData'',0);','FILENAME1973 = get(OH1973,''string'');','if isempty(FILENAME1973);','fprintf(''Filename Missing'');','else;','channames1973 = zeros(chans1973,6);','for i1973 = 1:chans1973;','Aa1973 = findobj(''Style'',''edit'',''UserData'',i1973);','elabel1973 = get(Aa1973,''string'');','channames1973(i1973,6-length(elabel1973)+1:6) = elabel1973;','end;','for i1973 = 1:length(channames1973(:));','if channames1973(i1973) == 0;','channames1973(i1973) = ''.'';','end;','end;','FID1973 = fopen(FILENAME1973,''w'');','fprintf(FID1973,''%s'',channames1973'');','fclose(FID1973);','fprintf(''Saving file'');','end;','clear chans1973 OH1973 FILENAME1973  Aa1973 elabel1973 ii1973 FID1973 channames1973 i1973;'];   
saveb = uicontrol('Units','Pixels','Style','PushButton','Position',[fl/10 20 fl*.15 25],'String','Save','UserData',chans,'Callback',TIMESTRING);

TIMESTRING = ['global is_done_pressed;','is_done_pressed = 1;','clear is_done_pressed'];
closeb = uicontrol('Units','Pixels','Style','PushButton','Position',...
               [fl*.45 20 fl*.15 25],'String','Done','Callback',TIMESTRING);

TIMESTRING = ['global is_done_pressed;','is_done_pressed = 2;','clear is_done_pressed'];
cancelb = uicontrol('Units','Pixels','Style','PushButton','Position',...
               [fl*.75 20 fl*.15 25],'String','Cancel','Callback',TIMESTRING);

while(is_done_pressed == 0)
   drawnow;
end

if is_done_pressed == 1
   channames = [];
   for i = 1:chans
      a = findobj('Style','edit','UserData',i);
      elabel = get(a,'string');
      channames = str2mat(channames,elabel);
   end
   channames = str2mat(channames,' ');   
else
 
   channames = [];
end

delete(fighandle)
clear is_done_pressed
