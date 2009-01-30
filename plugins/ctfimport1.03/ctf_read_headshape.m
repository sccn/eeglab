function HeadShape = ctf_read_headshape(file)

% ctf_read_headshape - read a CTF .shape file
%
% HeadShape = ctf_read_headshape(fileName)
%
% The *.shape file is an ascii text file in the following format:
%
% Number of Points
% x1 y1 z1
% x2 y2 z2
% .
% .
% .
% xn yn zn
%
%
% These vertex coordinates are returned into HeadShape (Nx3).  The
% coordinate values are in centimeters in either the voxel MRI
% coordinate system or the MEG Head Coordinate System (see
% ctf_read_mri for more about the coordinate system).
%

% $Revision: 1.1 $ $Date: 2009-01-30 03:49:27 $

% Copyright (C) 2004  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% History:  02/2004, Darren.Weber_at_radiology.ucsf.edu
%                    - adapted from an appendex to CTF document
%                    MRIViewer.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ver = '[$Revision: 1.1 $]';
fprintf('CTF_READ_HEADSHAPE [v%s]\n',ver(12:16));  tic;

fprintf('...checking file input parameter\n');

if ~exist('file','var'),
  [fileName, filePath, filterIndex] = uigetfile('*.shape', 'Locate CTF .shape file');
  file = fullfile(filePath, fileName);
elseif isempty(file),
  fprintf('...file is empty\n');
  [fileName, filePath, filterIndex] = uigetfile('*.shape', 'Locate CTF .shape file');
  file = fullfile(filePath, fileName);
end
if ~exist(file,'file'),
  fprintf('...file does not exist\n');
  [fileName, filePath, filterIndex] = uigetfile('*.shape', 'Locate CTF .shape file');
  file = fullfile(filePath, fileName);
end


fid = fopen(file,'r');

if isequal(fid,-1),
    S=sprintf('Could not open file: "%s"',file);
    error(S);
else
    
    fprintf('...reading CTF head shape file (.shape)\n');
    
    tic;
    
    % Check for comment on first line of file
    frewind(fid); temp = fscanf(fid,'%s',1); frewind(fid);
    if findstr(temp,'#'), temp = fgetl(fid); end
    
    % Read vertices
    Nvertices = fscanf(fid,'%d',1);
    fprintf('...reading %d vertices\n',Nvertices);
    vertices = fscanf(fid,'%f',[3,Nvertices]);
    % translate
    HeadShape = vertices';
    
    % Read faces
    %Nfaces = fscanf(fid,'%d',1);
    %fprintf('...Reading %d Faces\n',Nfaces);
    %faces = fscanf(fid,'%d',[4,Nfaces]);
    % remove first row (index) & translate
    %faces = faces(2:4,:)';
    
    fclose(fid);
    
    t = toc;
    fprintf('...done (%6.2f sec).\n\n',t);
    
end

return
