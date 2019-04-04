% floatwrite() - Write data matrix to float file.
%
% Usage: >> floatwrite(data,filename, 'format') 
%
% Inputs:
%   data     - write matrix data to specified file as four-byte floating point numbers.
%   filename - name of the file
%   'format' - The option FORMAT argument specifies the storage format as
%              defined by fopen. Default format is 'native'.
%   'transp|normal' - save the data transposed (.dat files) or not.
%
% Author: Sigurd Enghoff, CNL / Salk Institute, La Jolla, 7/1998 
%
% See also: floatread(), fopen()

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

% 07-08-99  FORMAT argument added -se
% 02-08-00  new version included in toolbox -sm
% 01-25-02 reformated help & license, added links -ad 

function A = floatwrite(A, fname, fform, transp)

if ~exist('fform')
	fform = 'native';
end
if nargin < 4
    transp = 'normal';
end

if strcmpi(transp,'normal')
    if strcmpi(class(A), 'mmo')
        A = changefile(A, fname);
        return;
    elseif strcmpi(class(A), 'memmapdata')
        % check file to overwrite
        % -----------------------
        [fpath1 fname1 ext1] = fileparts(fname);
        [fpath2 fname2 ext2] = fileparts(A.data.Filename);
        if isempty(fpath1), fpath1 = pwd; end
        
        fname1 = fullfile(fpath1, [fname1 ext1]);
        fname2 = fullfile(fpath2, [fname2 ext2]);
        if ~isempty(findstr(fname1, fname2))
            disp('Warning: raw data already saved in memory mapped file (no need to resave it)');
            return;
        end
        
        fid = fopen(fname,'wb',fform);
        if fid == -1, error('Cannot write output file, check permission and space'); end
        if size(A,3) > 1
            for ind = 1:size(A,3)
                tmpdata = A(:,:,ind);
                fwrite(fid,tmpdata,'float');
            end
        else
            blocks = [ 1:round(size(A,2)/10):size(A,2)];
            if blocks(end) ~= size(A,2), blocks = [blocks size(A,2)]; end
            for ind = 1:length(blocks)-1
                tmpdata = A(:, blocks(ind):blocks(ind+1));
                fwrite(fid,tmpdata,'float');
            end
        end
    else
        fid = fopen(fname,'wb',fform);
        if fid == -1, error('Cannot write output file, check permission and space'); end
        fwrite(fid,A,'float');
    end
else
    % save transposed
    for ind = 1:size(A,1)
        fwrite(fid,A(ind,:),'float');
    end
end;    
fclose(fid);
