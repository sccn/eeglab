function [p] = elec_write_emse(p)

% elec_write_emse - Write an EMSE probe file (*.elp)
% 
% Usage: [p] = elec_write_emse(p)
% 
% This script outputs x,y,z values to an EMSE 
% probe (*.elp) file.
% 
% EMSE *.elp files are in meters.  Also, when using
% EMSE *.elp files in the eeg_toolbox, it is necessary
% to swap X and Y.
% 
% See also: emse_read_elp, elec_open
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no express or implied warranties
% History:  10/2002, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[path,name,ext] = fileparts(strcat(p.elec.path,filesep,p.elec.file));
ext = '.elp';
file = fullfile(path,[name ext]);

p.elec.file = [name ext];

[fid,msg] = fopen(file,'w','ieee-le');
if ~isempty(msg), error(msg); end

tic

fprintf('\nELEC_WRITE_EMSE...\n');
fprintf('...writing .elp data to:\n\t%s\n',file);

write_elp(fid,p);

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [elp] = write_elp(fid,p)
    
    Nelec = size(p.elec.data.x,1);
    
    % Probe files contain position information for electrode locations 
    % and/or gradiometer locations. The file consists of a prolog, a 
    % header, and a list of one or more sensor fields.
    
    % Any line beginning with '//' is a comment line, which is ignored
    
    % Write the prolog
    fprintf(fid,'%d %d\n%d\n',3,2,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Write the header
    % The header consists of one optional entry and 2 entries in 
    % mandatory sequence and one optional entry:
    % Name [optional] > %N %s replace %s with name string (8 or fewer characters)
    % Type Code       > %x    replace %x with 1 (all electric), 2 (all magnetic) or 4 (mixed).
    % #Channels       > %d    number of points per channel per epoch [???? DLW]
    
    fprintf(fid,'//TypeCode nsensors\n');
    fprintf(fid,'%d %d\n',1,Nelec + 1);     % Nelectrodes plus REF
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fiducial points may be included optionally. They are required 
    % for MRI registration. If they are included, they must be in 
    % the obligatory order : nasion, left preauricular point, 
    % right preauricular point. Table A-2 defines the format for 
    % representing fiduciary points.
    
    fprintf(fid,'//Fiducials (Nasion, Left Preauricular, Right Preauricular)\n');
    fprintf(fid,'%s\t%12.8f\t%12.8f\t%12.8f\n','%F',...
        p.elec.data.nasion(1),...
        p.elec.data.nasion(2),...
        p.elec.data.nasion(3));
    fprintf(fid,'%s\t%12.8f\t%12.8f\t%12.8f\n','%F',...
        p.elec.data.lpa(1),...
        p.elec.data.lpa(2),...
        p.elec.data.lpa(3));
    fprintf(fid,'%s\t%12.8f\t%12.8f\t%12.8f\n','%F',...
        p.elec.data.rpa(1),...
        p.elec.data.rpa(2),...
        p.elec.data.rpa(3));
    
    % Each electrode is represented by an electric sensor, 
    % and consists of 5 fields, of which 1 (the name) is 
    % optional.
    % Name              Format      Description 
    % Type Code         %S          %x replace %x with 400 (electrode) or 1c00 if reference channel
    % Name [optional]   %N          %s replace %s with name string (8 or fewer characters)
    % Position          %g %g %g    electrode location with respect to head frame (Cartesian, meters)
    % Orientation       %g %g %g    not used, replace with 0 0 1
    % 
    % Sensor state (which appears in the 'type code' field) may 
    % be obtained by logically OR-ing suitable combinations from 
    % the table below. Note that not all combinations are physically valid.
    %
    % type/state        type code 
    % magnetic          200
    % electric          400
    % off               800
    % reference         1000        actually '1c00' [DLW]
    % optical           4000
    % trigger           8000
    % other             10000
    % named point       20000
    % 
    % Other types (such as named points, trigger, and optical) should 
    % be represented in the same pattern as electrodes, with the type 
    % code set to identify the type. Even those types (e.g. trigger) 
    % which do not have a true location, should have a nominal 
    % location, (e.g. 0 0 0).
    
    
    % Note below that eeg_toolbox x,y,z are equivalent to
    % emse y,-x,z
    
    
    typecode = 400;
    
    for n = 1:Nelec,
        fprintf(fid,'\n//ecSensor typecode/state---------------------\n');
        fprintf(fid,'%s\t%d\n','%S',typecode);
        fprintf(fid,'//ecSensor name:\n');
        fprintf(fid,'%s\t%s\n','%N',char(p.elec.data.label(n)));
        fprintf(fid,'//sphere origin\n');
        fprintf(fid,'%s\t%12.8f\t%12.8f\t%12.8f\n','%O',...
            p.elec.data.centroid(2),...
            p.elec.data.centroid(1) * -1,...
            p.elec.data.centroid(3));
        fprintf(fid,'//ecSensor location (origin)\n');
        fprintf(fid,'%s\t%12.8f\t%12.8f\t%12.8f\n','  ',...
            p.elec.data.y(n),...
            p.elec.data.x(n) * -1,...
            p.elec.data.z(n));
    end
    
    % Output the REF
    typecode = '1c00'; % This is not consistent with the
                       % table above, but it works!
    fprintf(fid,'\n//ecSensor typecode/state---------------------\n');
    fprintf(fid,'%s\t%s\n','%S',typecode);
    fprintf(fid,'//ecSensor name:\n');
    fprintf(fid,'%s\t%s\n','%N','Ref');
    fprintf(fid,'//sphere origin\n');
    fprintf(fid,'%s\t%12.8f\t%12.8f\t%12.8f\n','%O',...
        p.elec.data.centroid(2),...
        p.elec.data.centroid(1) * -1,...
        p.elec.data.centroid(3));
    fprintf(fid,'//ecSensor location (origin)\n');
    fprintf(fid,'%s\t%12.8f\t%12.8f\t%12.8f\n','  ',...
        p.elec.data.ref(2),...
        p.elec.data.ref(1) * -1,...
        p.elec.data.ref(3));
    
    fclose(fid);
    
return
