function elec_write_3dspace(p)

% elec_write_3dspace - Write electrode coordinate file (ascii).
% 
% Useage: elec_write_3dspace(p)
% 
% Write out the electrode labels, type, and Cartesian (x,y,z)
% coordinates in the format of Neuroscan 3Dspace ascii files.
% 
% Each row of the file comprises an electrode label, an 
% electrode type code (see below), and the x,y,z 
% coordinates (cm). Each field is separated by spaces.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no express or implied warranties
% History:  10/2002, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nELEC_WRITE_3DSPACE...\n'); tic;

if ~exist('p', 'var'),
    msg = sprintf('...no input p struct.\n\n');
    error(msg);
end

[path,name,ext] = fileparts(strcat(p.elec.path,filesep,p.elec.file));
ext = '.dat';
file = fullfile(path,[name ext]);

fprintf('...writing electrode data to:\n\t%s\n', file);

write_elec(file,p.elec.data);

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_elec(file,elec);

fid = fopen(file,'w','ieee-le');

if(fid == -1),
    fprintf('...could not open file:\n\t%s',file);
    return;
else
    
    % Write out the electrode coordinates
    
    fprintf('...converting from meters to cm.\n');
    % Convert to Neuroscan 3Dspace coordinate
    % metric (cm) from meters.
    
    nasion = elec.nasion .* 100;
    rpa    = elec.rpa .* 100;
    lpa    = elec.lpa .* 100;
    
    x = elec.x .* 100;
    y = elec.y .* 100;
    z = elec.z .* 100;
    
    ref = elec.ref .* 100;
    
    centroid = elec.centroid .* 100;
    
    
    fprintf('...writing output coordinates.\n');
    
    fprintf(fid,'%+10s\t%3d\t%+12.6f\t%+12.6f\t%+12.6f\n','Nasion',  110, nasion);
    fprintf(fid,'%+10s\t%3d\t%+12.6f\t%+12.6f\t%+12.6f\n','Left',    108, lpa);
    fprintf(fid,'%+10s\t%3d\t%+12.6f\t%+12.6f\t%+12.6f\n','Right',   114, rpa);
    
    Nelec = size(x,1);
    type = 69;
    for e = 1:Nelec,
        fprintf(fid,'%+10s\t%3d\t%+12.6f\t%+12.6f\t%+12.6f\n',char(elec.label(e)), type, x(e), y(e), z(e));
    end
    
    fprintf(fid,'%+10s\t%3d\t%+12.6f\t%+12.6f\t%+12.6f\n','Centroid', 99, centroid(1),centroid(2),centroid(3));
    fprintf(fid,'%+10s\t%3d\t%+12.6f\t%+12.6f\t%+12.6f\n','Ref',     120, ref(1),ref(2),ref(3));
    
    fclose(fid);
    
end

return
    
