function ctf_write_sam_targets(vertices,file)

% ctf_write_sam_targets - write a CTF SAM target file
%
% ctf_write_sam_targets(vertices,file)
%
% vertices - is Nx3 list of coordinates in CTF head coordinates (in cm).
% file - an output path/file
%
% The output file is an ascii text file in the following format:
%
% Number of Points
% x1 y1 z1
% x2 y2 z2
% .
% .
% .
% xn yn zn
%
% These vertex coordinates are contained in the input vertices
% matrix (Nx3).  The coordinate values must be in centimeters in
% the MEG Head Coordinate system (see ctf_read_mri for more about 
% the coordinate system).
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
fprintf('CTF_WRITE_SAM_TARGETS [v%s]\n',ver(12:16));



%--------------------------------------------------
% output an ascii file

fid = fopen(file,'w');

if isequal(fid,-1),
    S=sprintf('Could not open file: "%s"',file);
    error(S);
else
    
    fprintf('...writing to file:\n   %s\n',file);
    fprintf('...writing CTF SAM target file, assuming head coordinates\n');
    tic;
    
    % Write vertices
    Nvertices = size(vertices,1);
    fprintf(fid,'%d\n',Nvertices);
    
    for v = 1:Nvertices,
      fprintf(fid,'%6.3f %6.3f %6.3f\n',vertices(v,1),vertices(v,2),vertices(v,3));
    end
    
    fclose(fid);
    
    t = toc;
    fprintf('...done (%6.2f sec).\n\n',t);
    
end

return
