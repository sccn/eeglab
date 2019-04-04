% gettext() - This function prints a dialog box on screen and waits for 
%             the user to enter a string. There is a cancel button which 
%             returns a value of [].
% Usage:
%   >> out = gettext(label1,label2,...,label7);
%
% Author: Colin Humphries, CNL / Salk Institute, La Jolla, 1997

% Copyright (C) 1997 Colin Humphries, Salk Institute 
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

% 01-25-02 reformated help & license -ad 

function out = gettext(label1,label2,label3,label4,label5,label6,label7);

global is_text_entered
is_text_entered = 0;

for i = 1:nargin
   eval(['leng(i) = length(label',num2str(i),');'])
end
figurelength = 15*1.0*max(leng);

if figurelength < 240
   figurelength = 290;
end

figureheight = 80+nargin*18;
% Set figure

figure('Position',[302 396 figurelength figureheight],'color',[.6 .6 .6],'NumberTitle','off')
f1 = gcf;

ax = axes('Units','Pixels','Position',[0 0 figurelength figureheight],'Visible','off','XLim',[0 figurelength],'YLim',[0 figureheight]);

for i = 1:nargin

   eval(['text(figurelength/2,',num2str(figureheight-(i-1)*18-20),',label',num2str(i),',''color'',''k'',''FontSize'',16,''HorizontalAlignment'',''center'')'])

end
% Set up uicontrols

TIMESTRING = ['global is_text_entered;is_text_entered = 1;clear is_text_entered'];

u = uicontrol('Style','Edit','Units','Pixels','Position',[15 20 figurelength-95 25],'HorizontalAlignment','left','Callback',TIMESTRING);

TIMESTRING = ['global is_text_entered;is_text_entered = 2;clear is_text_entered'];
v = uicontrol('Style','Pushbutton','Units','Pixels','Position',[figurelength-70 20 55 25],'String','Cancel','Callback',TIMESTRING);


while(is_text_entered == 0)
   drawnow
end


if is_text_entered == 1
   out = get(u,'string');
else
   out = [];
end


clear is_text_entered
delete(f1)

