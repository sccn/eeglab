function [p] = elec_open(p)

% elec_open - opens electrode data for the eeg_toolbox
% 
% Usage: [p] = elec_open( p )
% 
% p is the eeg_toolbox struct (see eeg_toolbox_defaults).
% If p is omitted, the function uses the defaults.
% 
% This function requires the fields:
% 
% p.elec.path  - the directory location of the file to load
% p.elec.file  - the name of the file to load
% p.elec.type  - the file format type
% p.elec.n     - how many electrodes to read
% p.elec.plot  - boolean, 1 = plot, 2 = no plot
% 
% The return values are in p.elec.data.  The metric of the
% electrode coordinates returned is meters.
% 
% Recognised file format types are:
% 
% 'cartesian','spherical1','spherical2'
% 'scan3ddasc'
% 'scantri'
% 'brainstorm'
% 'emse','elp'
% 
% See also: elec_load, elec_load_scan_3ddasc, elec_load_scan_tri,
%           elec_load_brainstorm, emse_read_elp
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no express or implied warranties
% History:  02/2002 Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('p','var'),[p] = eeg_toolbox_defaults; end
if isempty(p),[p] = eeg_toolbox_defaults; end

eegversion = '$Revision: 1.1 $';
fprintf('ELEC_OPEN [v %s]\n',eegversion(11:15));

[path,name,ext] = fileparts(strcat(p.elec.path,filesep,p.elec.file));
file = fullfile(path,[name ext]);

% Get electrode dataset, depeding on type of coordinates

electype = lower(p.elec.type);

switch electype,
  
  case {'cartesian','spherical1','spherical2'},
    
    [elec,type,X,Y,Z,theta,phi,r] = elec_load(file,p.elec.type,0,0,0,p.elec.n);
    % Get electrode centroid
    index = find(type == 99); 
    xo = X(index); 
    yo = Y(index); 
    zo = Z(index);
    p.elec.data.centroid = [xo yo zo];
    
    % Select electrodes only
    index = find(type == 69);
    elec = elec(index);
    x = X(index);
    y = Y(index);
    z = Z(index);
    theta = theta(index);
    phi   = phi(index);
    r     = r(index);
    
    p.elec.data.label = elec;
    p.elec.data.x = x;
    p.elec.data.y = y;
    p.elec.data.z = z;
    p.elec.data.theta = theta;
    p.elec.data.phi   = phi;
    p.elec.data.r     = r;
    
    index = find(type == 110);
    p.elec.data.nasion = [X(index) Y(index) Z(index)];
    index = find(type == 108);
    p.elec.data.lpa    = [X(index) Y(index) Z(index)];
    index = find(type == 114);
    p.elec.data.rpa    = [X(index) Y(index) Z(index)];
    index = find(type == 120);
    p.elec.data.ref    = [X(index) Y(index) Z(index)];
    
    
    
  case {'emse','elp'},
    
    elp = emse_read_elp(file);
    
    % an example elp struct:
    %        version: 3
    %       filetype: 2
    %      minor_rev: 1
    %     sensorType: 4001
    %        sensorN: 125
    %         nasion: [0.0957 0 0]
    %            lpa: [-7.1503e-004 0.0804 0]
    %            rpa: [7.1503e-004 -0.0804 0]
    %              x: [124x1 double]
    %              y: [124x1 double]
    %              z: [124x1 double]
    %            ref: [0.0089 -0.0732 -0.0214]
    %         origin: [-0.0083 0.0043 0.0496]
    %           type: {124x1 cell}
    %           name: {124x1 cell}
    
    % EMSE coordinate orientation is +X anterior and +Y left,
    % whereas eeg_toolbox is         +Y anterior and +X right
    % effectively rotated -90 degrees
    
    p.elec.data.ref(1) = elp.ref(2) * -1;
    p.elec.data.ref(2) = elp.ref(1);
    p.elec.data.ref(3) = elp.ref(3);
    
    p.elec.data.x = elp.y * -1;
    p.elec.data.y = elp.x;
    p.elec.data.z = elp.z;
    
    p.elec.data.centroid(1) = elp.origin(2) * -1;
    p.elec.data.centroid(2) = elp.origin(1);
    p.elec.data.centroid(3) = elp.origin(3);
    
    p.elec.data.nasion  = elp.nasion;
    p.elec.data.lpa     = elp.lpa;
    p.elec.data.rpa     = elp.rpa;
    
    p.elec.data.label   = elp.name;
    
    
  case 'brainstorm',
    
   [p] = elec_load_brainstorm(p);
    
    
  case 'scan3ddasc',
    
    scan3dd = elec_load_scan_3ddasc(file);
    
    fprintf('...converting electrode XYZ from cm to meters.\n');
    
    p.elec.data.x = scan3dd.x ./ 100;
    p.elec.data.y = scan3dd.y ./ 100;
    p.elec.data.z = scan3dd.z ./ 100;
    
    p.elec.data.nasion  = scan3dd.nasion ./ 100;
    p.elec.data.lpa     = scan3dd.lpa ./ 100;
    p.elec.data.rpa     = scan3dd.rpa ./ 100;
    
    p.elec.data.label   = scan3dd.label;
    
    p.elec.data.ref = scan3dd.ref ./ 100;
    
    p.elec.data.centroid = scan3dd.origin ./ 100;
    
    
  case 'scantri',
    
    % Neuroscan 3Dspace TRI file
    
    tri = elec_load_scan_tri(file);
    
    p.elec.data.label = tri.label;
    
    fprintf('...converting electrode XYZ from cm to meters.\n');
    p.elec.data.x = tri.XYZ(:,1) ./ 100;
    p.elec.data.y = tri.XYZ(:,2) ./ 100;
    p.elec.data.z = tri.XYZ(:,3) ./ 100;
    
    fprintf('...guessing origin is at (0,0,0).\n');
    p.elec.data.centroid = [0 0 0];
    
    
  otherwise,
    
    msg = sprintf('...cannot read file types %s\n',electype);
    error(msg);
    
end



% -- Calculate some extra parameters

x = p.elec.data.x;
y = p.elec.data.y;
z = p.elec.data.z;

xo = p.elec.data.centroid(1);
yo = p.elec.data.centroid(2);
zo = p.elec.data.centroid(3);

if ~isfield(p.elec.data,'theta'),
  [theta,phi,r] = elec_cart2sph(x,y,z,xo,yo,zo);	
  p.elec.data.theta = theta;
  p.elec.data.phi   = phi;
  p.elec.data.r     = r;
end

% Estimate X,Y,Z radii
Xrad = (max(x)-min(x))/2; 
Yrad = (max(y)-min(y))/2; 
Zrad = (max(z)-min(z));
p.elec.data.R = [Xrad Yrad Zrad];

% Estimate ellipse that best fits electrode co-ordinates
[r,Xel,Yel,Zel] = elec_ellipse_fit(x,y,z,xo,yo,zo,100,p.elec.plot);
p.elec.data.Xel = Xel;
p.elec.data.Yel = Yel;
p.elec.data.Zel = Zel;
p.elec.data.Rel = r;

% Estimate sphere that best fits electrode co-ordinates
[r,Xsp,Ysp,Zsp] = elec_sphere_project(x,y,z,xo,yo,zo,1,p.elec.plot);
p.elec.data.Xsp = Xsp;
p.elec.data.Ysp = Ysp;
p.elec.data.Zsp = Zsp;
p.elec.data.Rsp = [r r r];

% Create a refined spherical mesh and the interpolation
% matrices - used in topographic mapping
%p = elec_sph_refine(p);
%p.elec.data.Lsp = mesh_laplacian(p.elec.data.Vsp,p.elec.data.Fsp);
%p.elec.data.Isp = mesh_laplacian_interp(p.elec.data.Lsp, 1:length(p.elec.data.Xsp));

% Define the electrode regions
p.elec.data.regions = elec_regions;

return
