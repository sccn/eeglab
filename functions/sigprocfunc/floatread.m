% floatread() - Read matrix from float file ssuming four byte floating point number
%
% Usage: >> A = floatread(filename,size,'format',offset) 
%
% Inputs:
%   filename - name of the file
%   size     - determine the number of float elements to be read and 
%              the dimensions of the resulting matrix. If the last element 
%              of 'size' is INF the size of the last dimension is determined
%              by the file length. If size is 'square,' floatread() attempts 
%              to read a square matrix.
%
% Optional inputs:
%  'format'  - the option FORMAT argument specifies the storage format as
%              defined by fopen(). Default format ([]) is 'native'.
%  offset    - offset in floats from the beginning of file (=0)
%              to start reading (4-byte floats assumed).
%              It uses fseek to read a portion of a large data file.
%
% Author: Sigurd Enghoff, CNL / Salk Institute, La Jolla, 7/1998 
%
% See also: floatwrite(), fopen()

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

% 04-26-99  modified by Sigurd Enghoff to handle variable-sized and
%           multi-dimensional data.
% 07-08-99  modified by Sigurd Enghoff, FORMAT argument added.
% 02-08-00  help updated for toolbox inclusion -sm
% 02-14-00  added segment arg -sm
% 08-14-00  added size 'square' option -sm
% 01-25-02  reformated help & license, added links -ad 

function A = floatread(fname,Asize,fform,offset)

if nargin<2
  help floatread
  return
end

if ~exist('fform') | isempty(fform)|fform==0
	fform = 'native';
end

fid = fopen(fname,'rb',fform);
if fid>0 
 if exist('offset')
   stts = fseek(fid,4*offset,'bof');
   if stts ~= 0
     error('floatread(): fseek() error.');
     return
   end
 end
 if ischar('Asize')
   if strcmp(Asize,'square')
         fseek(fid,0,'eof'); % go to end of file
         bytes = ftell(fid); % get byte position
         fseek(fid,0,'bof'); % rewind
         bytes = bytes/4; % nfloats
         froot = sqrt(bytes);
         if round(froot)*round(froot) ~= bytes
              error('floatread(): filelength is not square.')
         else
              Asize = [round(froot) round(froot)];
         end
   end
 end
 A = fread(fid,prod(Asize),'float');
else
 error('floatread() fopen() error.');
 return
end
% fprintf('   %d floats read\n',prod(size(A)));

if Asize(end) == Inf
	Asize = Asize(1:end-1);
	A = reshape(A,[Asize length(A)/prod(Asize)]);
else
	A = reshape(A,Asize);
end

fclose(fid);
