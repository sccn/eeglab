% shortread() - Read matrix from short file.
%
% Usage:
%  >> A = shortread(filename,size,'format',offset) 
%
% Inputs:
%  filename - Read matrix a from specified file while assuming four byte  
%             short integers. 
%  size     - The vector SIZE determine the number of short elements to be 
%             read and the dimensions of the resulting matrix. If the last 
%             element of SIZE is INF the size of the last dimension is determined 
%             by the file length.
%
% Optional inputs:
% 'format'  - The option FORMAT argument specifies the storage format as
%             defined by fopen(). Default format ([]) is 'native'.
% offset    - The option OFFSET is offset in shorts from the beginning of file (=0)
%             to start reading (2-byte shorts assumed).
%             It uses fseek to read a portion of a large data file.
%
% Author: Sigurd Enghoff, CNL / Salk Institute, La Jolla, 7/1998 
%
% See also: floatread(), floatwrite(), fopen()

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

% xx-07-98  written by Sigurd Enghoff, Salk Institute
% 04-26-99  modified by Sigurd Enghoff to handle variable-sized and
%           multi-dimensional data.
% 07-08-99  modified by Sigurd Enghoff, FORMAT argument added.
% 02-08-00  help updated for toolbox inclusion -sm
% 02-14-00  added segment arg -sm

function A = shortread(fname,Asize,fform,offset)

if ~exist('fform') || isempty(fform) || fform(1)==0
	fform = 'native';
end

fid = fopen(fname,'rb',fform);
if fid>0 
 if exist('offset')
   stts = fseek(fid,2*offset,'bof');
   if stts ~= 0
     error('shortread(): fseek() error.');
     return
   end
 end
 A = fread(fid,prod(Asize),'short');
else
 error('shortread(): fopen() error.');
 return
end
fprintf('   %d shorts read\n',prod(size(A)));

if Asize(end) == Inf
	Asize = Asize(1:end-1);
	A = reshape(A,[Asize length(A)/prod(Asize)]);
else
	A = reshape(A,Asize);
end

fclose(fid);
