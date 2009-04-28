function [scan3dd] = elec_load_scan_3ddasc(filename)

% elec_load_scan_3ddasc - read ascii export of Neuroscan 3DD file
% 
% Load an ascii electrode file and return 
% the electrode labels and coordinates.  This
% function can read and return Cartesian (x,y,z) 
% and/or spherical (theta,phi,r) coordinates.
%
% [scan3dd] = elec_load_scan_3ddasc(filename)
% 
% where:
% 
% file = '<path><filename>' with row format '%s %d %f %f %f'
% 
% The file format is that of NeuroScan 3Dspace ascii 
% export files.  Each row of the file comprises an electrode label, 
% an electrode type code (see below), and the x,y,z coordinates (cm).
% Each field is separated by spaces.  For example,
% 
% Fz    69  Xcm Ycm Zcm
% 
% Example result:
%
% scan3dd = 
% 
%      label: {1x128 cell}
%          x: [128x1 double]
%          y: [128x1 double]
%          z: [128x1 double]
%        hsp: [1617x3 double]
%     nasion: [-0.0333 9.0341 0]
%        lpa: [-6.9128 0 0]
%        rpa: [6.9128 0 0]
%     origin: [0 0 0]
%        ref: [0 9.9341 0]
%
% All coordinates are in centimeters; +X is right, +Y is anterior
% and +Z is superior.
% 
% Notes:
% 
% i) Type is defined as follows (from Neuroscan 3Dspace ascii export):
%
%           Electrode               Type
%           ---------               ----
%           Nasion                   110 (or 78)
%           Left                     108 (or 76)
%           Right                    114 (or 82)
%           electrodes                69
%           Centroid (origin)         99 (or 67; eg <0,0,0>)
%           Ref                      120 (or 88)
%           Scalp Points (hsp)        32
%
% 


% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no express or implied warranties
% History:  08/1999, Darren.Weber_at_radiology.ucsf.edu
%           08/2003, Darren.Weber_at_radiology.ucsf.edu
%                    complete rewrite, replacing elec_load.m
%                    new version of 3Dspace exports different type numbers ;-)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ver = '$Revision: 1.1 $';
fprintf('\nELEC_LOAD_SCAN_3DDASC [v %s]\n',ver(11:15));

tic;

[path,name,ext] = fileparts(filename);
file = fullfile(path,[name ext]);

fprintf('...loading electrodes from:\n\t%s\n', file);

scan3dd = read_3dd(file);

fprintf('...loaded %d electrodes\n', size(scan3dd.x,1));

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scan3dd = read_3dd(file),


scan3dd.label = [];
scan3dd.x = [];
scan3dd.y = [];
scan3dd.z = [];
scan3dd.hsp = [];



fid = fopen(file);
if fid < 0,
  msg = sprintf('cannot open file: %s\n',file);
  error(msg);
end

% 3DD ascii files contain position information for 
% fiducial, sensor, and head shape fields.

% First get the fiducial points by reading the whole file
% (clumsy, but exhaustive search for field indicators).

% Fiducial points are required for MRI registration. They are 
% the nasion, left and right preauricular points, eg:

%      Nasion	78	-0.033340	9.034104	0.000000
%      Left	76	-6.912818	-0.000000	0.000000
%      Right	82	6.912818	0.000000	-0.000000

fprintf('...searching for fiducials...');

n = 0;
while n < 5,
  tmp = fgetl(fid);
  if tmp < 0, break; end
  
  tmp = lower(tmp);
  
  if strfind(tmp,'nasion'),
    tmp = sscanf(tmp,'%s %d %f %f %f');
    scan3dd.nasion = [tmp(end-2) tmp(end-1) tmp(end)];
    n = n + 1;
    continue;
  end
  if strfind(tmp,'left'),
    tmp = sscanf(tmp,'%s %d %f %f %f');
    scan3dd.lpa = [tmp(end-2) tmp(end-1) tmp(end)];
    n = n + 1;
    continue;
  end
  if strfind(tmp,'right'),
    tmp = sscanf(tmp,'%s %d %f %f %f');
    scan3dd.rpa = [tmp(end-2) tmp(end-1) tmp(end)];
    n = n + 1;
    continue;
  end
  if strfind(tmp,'centroid'),
    tmp = sscanf(tmp,'%s %d %f %f %f');
    scan3dd.origin = [tmp(end-2) tmp(end-1) tmp(end)];
    n = n + 1;
    continue;
  end
  if strfind(tmp,'ref'),
    tmp = sscanf(tmp,'%s %d %f %f %f');
    scan3dd.ref = [tmp(end-2) tmp(end-1) tmp(end)];
    n = n + 1;
    continue;
  end
end

frewind(fid);
fprintf('done\n');


fprintf('...searching for electrodes...');

ok = 1;
while ok,
  tmp = fgetl(fid);
  if tmp < 0, break; end
  
  if strfind(lower(tmp),'nasion'),   continue; end
  if strfind(lower(tmp),'left'),     continue; end
  if strfind(lower(tmp),'right'),    continue; end
  if strfind(lower(tmp),'centroid'), continue; end
  if strfind(lower(tmp),'ref'),      continue; end
  
  tmp = sscanf(tmp,'%s %d %f %f %f');
  
  if tmp(end-3) == 69,
    scan3dd.label{end+1} = char(tmp(1:end-4))';
    scan3dd.x(end+1,1) = tmp(end-2);
    scan3dd.y(end+1,1) = tmp(end-1);
    scan3dd.z(end+1,1) = tmp(end);
    continue;
  end
  if strfind(char(tmp(1:end-4))','32'),
    break;
    % found a head shape point
    scan3dd.hsp(end+1,:) = [tmp(end-2),tmp(end-1),tmp(end)];
    continue;
  end
end


frewind(fid);
fprintf('done\n');



fprintf('...searching for head shape points...');

ok = 1;
while ok,
  tmp = fgetl(fid);
  if tmp < 0, break; end
  
  if strfind(lower(tmp),'nasion'),   continue; end
  if strfind(lower(tmp),'left'),     continue; end
  if strfind(lower(tmp),'right'),    continue; end
  if strfind(lower(tmp),'centroid'), continue; end
  if strfind(lower(tmp),'ref'),      continue; end
  if strfind(lower(tmp),'69'),       continue; end
  
  tmp = sscanf(tmp,'%d %f %f %f');
  
  if tmp(end-3) == 32,
    % found a head shape point
    scan3dd.hsp(end+1,:) = [tmp(end-2),tmp(end-1),tmp(end)];
    continue;
  end
end

frewind(fid);
fprintf('done\n');

fclose(fid);


return
