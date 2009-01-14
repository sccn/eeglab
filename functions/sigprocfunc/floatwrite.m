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
% Revision 1.8  2008/11/24 21:51:17  arno
% prevent overwriting memory map file
%
% Revision 1.7  2008/11/22 02:59:00  arno
% use temporary variable
%
% Revision 1.6  2008/11/19 23:26:02  arno
% save as transposed if necessary
%
% Revision 1.5  2008/11/17 19:34:55  ywu
% same
%
% Revision 1.4  2008/11/17 19:33:46  ywu
% same
%
% Revision 1.3  2008/11/17 19:30:07  ywu
% same
%
% Revision 1.2  2008/11/17 19:26:34  ywu
% support for memory map data
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 07-08-99  FORMAT argument added -se
% 02-08-00  new version included in toolbox -sm
% 01-25-02 reformated help & license, added links -ad 

function floatwrite(A, fname, fform, transp)

if ~exist('fform')
	fform = 'native';
end
if nargin < 4
    transp = 'normal';
end;

if strcmpi(transp,'normal')
    if strcmpi(class(A), 'memmapdata')
        
        % check file to overwrite
        % -----------------------
        [fpath1 fname1 ext1] = fileparts(fname);
        [fpath2 fname2 ext2] = fileparts(A.data.Filename);
        if isempty(fpath1), fpath1 = pwd; end;
        
        fname1 = fullfile(fpath1, [fname1 ext1]);
        fname2 = fullfile(fpath2, [fname2 ext2]);
        if ~isempty(findstr(fname1, fname2))
            disp('Warning: data already saved (cannot overwrite memory map file)');
            return;
        end;
        
        fid = fopen(fname,'wb',fform);
        if fid == -1, error('Cannot write output file, check permission and space'); end;
        if size(A,3) > 1
            for ind = 1:size(A,3)
                tmpdata = A(:,:,ind);
                fwrite(fid,tmpdata,'float');
            end;
        else
            blocks = [ 1:round(size(A,2)/10):size(A,2)];
            if blocks(end) ~= size(A,2), blocks = [blocks size(A,2)]; end;
            for ind = 1:length(blocks)-1
                tmpdata = A(:, blocks(ind):blocks(ind+1));
                fwrite(fid,tmpdata,'float');
            end;
        end;
    else
        fid = fopen(fname,'wb',fform);
        if fid == -1, error('Cannot write output file, check permission and space'); end;
        fwrite(fid,A,'float');
    end;
else
    % save transposed
    for ind = 1:size(A,1)
        fwrite(fid,A(ind,:),'float');
    end;
end;    
fclose(fid);
