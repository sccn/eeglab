function [tri] = elec_load_scan_tri(file)

% elec_load_scan_tri - read sensor coordinates from Scan .tri files
% 
% INPUTS:
%     file  - string, name of the .tri file
% OUTPUTS:
%     tri.hdr   - header and electrode fields
%     tri.XYZ   - electrode coordinates (tri.hdr.elec)
%     tri.label - first 4 letters of electrode names
% 
% Also check small script at the end of the code to print
% out coordinates together with names (a fast hack)
% 
% Information on Neuroscan .tri format is at
% http://www.neuro.com/neuroscan/triformat.htm, which
% is copied at the end of this function .m file
% 
% .tri faces and vertices are not returned in this version
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence: GNU GPL, no express or implied warranties
% History: 10/2002, Yaroslav Halchenko, CS Dept. UNM
%                   (yoh@onerussian.com, ICQ#: 60653192)
%          01/2003, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(file,'r');

if isequal(fid,-1),
    msg = sprintf('Could not open file: "%s"',file);
    error(msg);
else
    
    fprintf('...Reading NeuroScan Tesselation (.tri)');
    
    tri.hdr.ID       = fread(fid,1,'long');     % Long ID (should be 100003, or 100004. (or 100002) )
    tri.hdr.filetype = fread(fid,1,'short');    % Short Filetype (=2 for triangle file)
    tri.hdr.rev      = fread(fid,1,'short');    % short revision
    tri.hdr.elthick  = fread(fid,1,'float');    % float electrodethickness
    tri.hdr.eldiam   = fread(fid,1,'float');    % float electrodediameter
    tri.hdr.reserved = fread(fid, 4080,'char'); % BYTE reserved[4080]
    
    % Then the face and vertex information follows:
    
    tri.hdr.Nfaces     = fread(fid,1,'short');
    tri.hdr.Nvertices  = fread(fid,1,'short');
    
    % Then for all the faces, 
    % the centroid (centre of face) x,y,z coordinates (unit vector)
    % and it's length are writen as four floats
    centroid = fread(fid,4* tri.hdr.Nfaces,'float');
%     tri.centroid = reshape(centroid,tri.hdr.Nfaces,4);
    
    % Then follows the vertex coordinates,
    % x, y, z, (normalized) and it's length
    vertices = fread(fid,4* tri.hdr.Nvertices,'float');
    
%     vertices = reshape(vertices,tri.hdr.Nvertices,4);
%     tri.vert = zeros(tri.hdr.Nvertices,3);
%     tri.vert(:,1) = vertices(:,1) .* vertices(:,4);
%     tri.vert(:,2) = vertices(:,2) .* vertices(:,4);
%     tri.vert(:,3) = vertices(:,3) .* vertices(:,4);
    
    % Then for all the faces, the three vertices that belong to them
    faces = fread(fid,3* tri.hdr.Nfaces,'short');
%     tri.faces = reshape(faces,tri.hdr.Nfaces,3) + 1;  % add 1 for matlab
    
    %Then the number of electrodes follows
    tri.hdr.Nelectrodes = fread(fid,1,'ushort');
    
    % Then for all electrodes
    for e = 1:tri.hdr.Nelectrodes,
        tri.hdr.elec(e).label = fread(fid,10,'char')'; % label of electrode, max 9 chars + \0
        tri.hdr.elec(e).key   = fread(fid,1,'short');  % key, normally = 'e' for electrode
        tri.hdr.elec(e).pos   = fread(fid,3,'float')'; % x, y, z (position)
        tri.hdr.elec(e).index = fread(fid,1,'ushort'); % electrode index number
    end
    
    fclose(fid);
    
    tmp = char(tri.hdr.elec(:).label);
    tri.label = tmp(:,1:4);
    
    tri.XYZ = reshape([tri.hdr.elec(:).pos],3,tri.hdr.Nelectrodes)';
    
%     tri.XYZ = zeros(e,3);
%     for e = 1:tri.hdr.Nelectrodes
%         tri.XYZ(e,1:3) = tri.hdr.elec(e).pos;
%     end
    
end

return;


%----------------------------------------------------------------

fout = fopen('XYZ.locs','w');
for e = 1:tri.hdr.Nelectrodes,
    nop = 0;
    for z=1:10,
        if (tri.hdr.elec(e).label(z)==0), nop=1; end
        if (~nop),
            fprintf(fout,'%c',tri.hdr.elec(e).label(z));
        else
            fprintf(fout,' ');
        end
    end
    
    fprintf(fout, '%8.3f ',  tri.hdr.elec(e).pos);
    fprintf(fout, '\n');
end
fclose(fout);

return

% --------------------------------------------------------------------------------


%  3D Space TRI file format
% 
% Below you will find the description of the TRI 
% file format as used in the 3D Space program. The 
% TRI file contains the description of a triangulated 
% head surface. The description uses C/C++ variables 
% to describe elements.
% 
% 
% --------------------------------------------------------
% 
% Each file starts with a header, which looks like this:
% 
% Long ID (should be 100003, or 100004. (or 100002) )
% Short Filetype (=2 for triangle file)
% short revision
% float electrodethickness
% float electrodediameter
% BYTE reserved[4080]
% 
% Then the facet and vertex information follows:
% 
% short number_of_facets
% short number_of_vertices
% 
% Then for all the facets give by number_of_facets:
% 
% the centroid (centre of facet) x,y,z coordinates
% (unit vector) and it's length are writen as four floats. 
% 
% Then follows the facet vertex coordinates for all 
% vertices (given by number_of_vertices):
% 
% four floats: x, y, z, (normalized) and it's length
% 
% Then for all the facets (give by number_of_facets):
% 
% the three vertices that belong to this facet. These 
% are three shorts, and are index values to the proper 
% vertice as listed above. (largest number is given by 
% number_of_vertices).
% 
% Then the number of electrodes follows:
% 
% (unsigned short) number_of_electrodes
% 
% Then for all electrodes (given by number_of_electrodes):
% 
% char label[10] (label of electrode, max 9 chars + \0)
% short key (normally = 'e' for electrode)
% float x, y, z (position)
% unsigned short ix (electrode index number)
% 
