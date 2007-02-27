% errordlg2() - Makes a popup dialog box with the specified message and (optional)
%               title.
%
% Usage:
%   errordlg2(Prompt, Title);
%
% Example:
%   errordlg2('Explanation of error','title of error');
%
% Input:
%   Prompt  -   A text string explaning why the user is seeing this error message.
%   Title   _   A text string that appears in the title bar of the error message.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 12 August 2002
%
% See also: inputdlg2(), questdlg2()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, arno@salk.edu
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

% $Log: not supported by cvs2svn $
% Revision 1.4  2002/11/15 02:15:07  arno
% header typos
%
% Revision 1.3  2002/08/28 01:04:34  arno
% debugging beep
%
% Revision 1.2  2002/08/14 16:44:54  arno
% beep update
%
% Revision 1.1  2002/08/12 18:48:22  arno
% Initial revision
%
% Revision 1.3  2002/08/12 18:24:29  arno
% debug
%
% Revision 1.2  2002/08/12 18:02:47  arno
% debug
%
% Revision 1.1  2002/08/12 18:01:34  arno
% Initial revision
%

function errordlg2(Prompt, Title);

if exist('beep') == 5
	beep;
else
	disp(char(7));
end;
if nargin <2
	Title = 'Error';
end;
questdlg2(Prompt, Title, 'OK', 'OK');