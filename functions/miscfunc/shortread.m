%shortread() - Read matrix from short file.
%
%   >> A = shortread(filename,size,[format],[offset]) 
%
%   Read matrix a from specified
%   file while assuming four byte short integers. The vector
%   SIZE determine the number of short elements to be read and the dimensions
%   of the resulting matrix. If the last element of SIZE is INF the
%   size of the last dimension is determined by the file length.
%
%   The option FORMAT argument specifies the storage format as
%   defined by fopen(). Default format ([]) is 'native'.
%
%   The option OFFSET is offset in shorts from the beginning of file (=0)
%   to start reading (2-byte shorts assumed).
%   It uses fseek to read a portion of a large data file.
%
%   See also floatread(), floatwrite(), fopen.

% xx-07-98  written by Sigurd Enghoff, Salk Institute
% 04-26-99  modified by Sigurd Enghoff to handle variable-sized and
%           multi-dimensional data.
% 07-08-99  modified by Sigurd Enghoff, FORMAT argument added.
% 02-08-00  help updated for toolbox inclusion -sm
% 02-14-00  added segment arg -sm

function A = shortread(fname,Asize,fform,offset)

if ~exist('fform') | isempty(fform)|fform==0
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
