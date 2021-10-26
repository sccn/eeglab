% floatread() - Read matrix from float file ssuming four byte floating point number
%               Can use fseek() to read an arbitrary (contiguous) submatrix.
%
% Usage:        >> a = floatread(filename,size,'format',offset) 
%
% Inputs:
%   filename - name of the file
%   size     - determine the number of float elements to be read and 
%              the dimensions of the resulting matrix. If the last element 
%              of 'size' is Inf, the size of the last dimension is determined
%              by the file length. If size is 'square,' floatread() attempts 
%              to read a square 2-D matrix.
%
% Optional inputs:
%  'format'  - the option FORMAT argument specifies the storage format as
%              defined by fopen(). Default format ([]) is 'native'.
%  offset    - either the number of first floats to skip from the beginning of the
%              float file, OR a cell array containing the dimensions of the original 
%              data matrix and a starting position vector in that data matrix. 
%
%              Example: % Read a [3 10] submatrix of a four-dimensional float matrix 
%                >> a = floatread('mydata.fdt',[3 10],'native',{[[3 10 4 5],[1,1,3,4]});
%              % Note: The 'size' and 'offset' arguments must be compatible both
%              % with each other and with the size and ordering of the float file.
% 
% Author: Sigurd Enghoff, CNL / Salk Institute, La Jolla, 7/1998 
%
% See also: floatwrite(), fopen()

% Copyright (C) Sigurd Enghoff, CNL / Salk Institute, La Jolla, 7/1998
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

if ~exist('fform') || isempty(fform) || fform(1)==0
	fform = 'native';
end

if ~exist('offset') 
	offset = 0;
end

fid = fopen(fname,'rb',fform);
if fid>0 
 if exist('offset')
   if iscell(offset)
     if length(offset) ~= 2
        error('offset must be a positive integer or a 2-item cell array');
     end
     datasize = offset{1};
     startpos = offset{2};
     if length(datasize) ~= length(startpos)
        error('offset must be a positive integer or a 2-item cell array');
     end
     for k=1:length(datasize)
       if startpos(k) < 1 || startpos(k) > datasize(k)
          error('offset must be a positive integer or a 2-item cell array');
       end
     end
     if length(Asize)> length(datasize)
        error('offset must be a positive integer or a 2-item cell array');
     end
     for k=1:length(Asize)-1
         if startpos(k) ~= 1 
            error('offset must be a positive integer or a 2-item cell array');
         end
     end
     sizedim = length(Asize);
     if Asize(sizedim) + startpos(sizedim) - 1 > datasize(sizedim)
        error('offset must be a positive integer or a 2-item cell array');
     end
     for k=1:length(Asize)-1
         if Asize(k) ~= datasize(k)
            error('offset must be a positive integer or a 2-item cell array');
         end
     end

     offset = 0;
     jumpfac = 1;
     for k=1:length(startpos)
           offset = offset + jumpfac * (startpos(k)-1);
           jumpfac = jumpfac * datasize(k);
     end

   elseif length(offset) > 1
     error('offset must be a positive integer or a 2-item cell array');
   end
   
   % perform the fseek() operation
   % -----------------------------
   stts = fseek(fid,4*offset,'bof');

   if stts ~= 0
     error('floatread(): fseek() error.');
     return
   end
 end

% determine what 'square' means
% -----------------------------
 if ischar('Asize')
   if iscell(offset)
         if length(datasize) ~= 2 || datasize(1) ~= datasize(2)
              error('size ''square'' must refer to a square 2-D matrix');
         end
         Asize = [datsize(1) datasize(2)];
   elseif strcmp(Asize,'square')
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

% interpret last element of Asize if 'Inf'
% ----------------------------------------
if Asize(end) == Inf
	Asize = Asize(1:end-1);
	A = reshape(A,[Asize length(A)/prod(Asize)]);
else
	A = reshape(A,Asize);
end

fclose(fid);
