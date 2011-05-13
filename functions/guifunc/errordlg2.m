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

function errordlg2(Prompt, Title);

if exist('beep') == 5
	beep;
else
	disp(char(7));
end;
if nargin <2
	Title = 'Error';
end;
if ~ismatlab, error(Prompt); end;
questdlg2(Prompt, Title, 'OK', 'OK');
