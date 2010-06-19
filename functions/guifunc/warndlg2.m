% warndlg2() - same as warndlg for eeglab()
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

function warndlg2(Prompt, Title);

if nargin <2
	Title = 'Warning';
end;
questdlg2(Prompt, Title, 'OK', 'OK');
