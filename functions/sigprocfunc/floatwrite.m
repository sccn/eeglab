% floatwrite() - Write data matrix to float file.
%
% Usage: >> floatwrite(data,filename, 'format') 
%
% Inputs:
%   data     - write matrix data to specified file as four-byte floating point numbers.
%   filename - name of the file
%   'format' - The option FORMAT argument specifies the storage format as
%              defined by fopen. Default format is 'native'.
%
% Author: Sigurd Enghoff, CNL / Salk Institute, La Jolla, 7/1998 
%
% See also: floatread(), fopen()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Sigurd Enghoff, CNL / Salk Institute, La Jolla, 7/1998
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

% 07-08-99  FORMAT argument added -se
% 02-08-00  new version included in toolbox -sm
% 01-25-02 reformated help & license, added links -ad 

function floatwrite(A,fname,fform)

if ~exist('fform')
	fform = 'native';
end

fid = fopen(fname,'wb',fform);
fwrite(fid,A,'float');
fclose(fid);
