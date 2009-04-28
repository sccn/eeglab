function [p] = mesh_scalp_interp(p)

% mesh_scalp_interp - Interpolates scalar over scalp mesh from electrodes
% 
% Usage: [p] = mesh_scalp_interp(p)
% 
% p is the eeg_toolbox struct (see eeg_toolbox_defaults)
% 
% This function returns the laplacian matrix of the scalp
% mesh, the interpolation matrix based on the min norm of
% the laplacian, and the voltage interpolated on the scalp.
% It also calls an electrode "coregistration" function and
% the results of that operation.
% 
% The return matrices are all in p.mesh.data, including the 
% interpolated timeseries values in p.mesh.data.timeseries{N},
% where N is the index of the scalp mesh in this case (see 
% mesh_check to get N).
% 
% This function is called by eeg_contours_engine.
% 
% See also, mesh_fit_elec, mesh_laplacian, mesh_laplacian_interp
% for more information.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  Gnu GPL, no express or implied warranties
% History:  09/02, Darren.Weber_at_radiology.ucsf.edu
%                  extracted out of MESH_SCALP_INTERP, called by that
%                  function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('p','var'),
    fprintf('Setting default parameters.\n');
   [p] = eeg_toolbox_defaults('create');
end

% Calculate voltage at scalp mesh vertices
if isempty(p.mesh.data),
    if ~isempty(p.mesh.file),
       [p] = mesh_open(p);
    else
        msg = sprintf('MESH_SCALP_INTERP: No meshes, load meshes first.');
        warning(msg);
    end
end

% get the scalp mesh
[p.mesh.current,meshExists] = mesh_check(p,'scalp');

if isempty(p.mesh.current),
    msg = sprintf('MESH_SCALP_INTERP: No ''scalp'' mesh, check p.mesh.data.meshtype.');
    error(msg);
end

% Find nearest vertices of scalp mesh to electrodes
% This will also reorder the scalp mesh vertices/faces
% so that the first 1:Nelec vertices are the vertices
% in the scalp that lie closest to the electrodes.
% Serious problems arise if there are any duplicate
% scalp vertices for any electrode, but this routine
% will continue regardless, at present.
p = mesh_fit_elec(p);

% Ensure scalp is the currently selected mesh
elecN = mesh_check(p,'elec');
scalpN = mesh_check(p,'scalp');
[p.mesh.current,meshExists] = mesh_check(p,'scalp');

% Get vertices and faces of scalp mesh, after they
% have been rearraned by mesh_fit_elec
scalpvert = p.mesh.data.vertices{scalpN};
scalpface = p.mesh.data.faces{scalpN};

% Get the scalp vertices and faces
% that correspond to the electrodes
elecvert  = p.mesh.data.vertices{elecN};
elecface  = p.mesh.data.faces{elecN};
% Get the index of each electrode vertex
% into the scalp mesh, should be first 
% 1:Nelec indices of scalp vertices
elecindex = p.mesh.data.elecindex{scalpN};

fprintf('MESH_SCALP_INTERP: Calculating Scalp Interpolation\n');


% Could identify the elements of the scalp below
% the electrodes and discard them before calculating
% the interpolation matrices, or use them in the 
% interpolation, but set their potential to zero
% later (see comment below).


% Calculate the Laplacian matrix of the scalp mesh (once off calculation)
if isfield(p.mesh.data,'lapmat'),
    if isempty(p.mesh.data.lapmat{scalpN}),
        p.mesh.data.lapmat{scalpN} = mesh_laplacian(scalpvert,scalpface);
    end
else
    p.mesh.data.lapmat{scalpN} = mesh_laplacian(scalpvert,scalpface);
end


% Calculate the electrode to scalp mesh interpolation matrix
% using a Laplacian minimum norm constraint (once off calculation)
if isfield(p.mesh.data,'lapint'),
    if isempty(p.mesh.data.lapint{scalpN}),
        [p.mesh.data.lapint{scalpN}, p.mesh.data.keepindex{scalpN}, p.mesh.data.repindex{scalpN}] = mesh_laplacian_interp(p.mesh.data.lapmat{scalpN}, elecindex);
    end
else
    [p.mesh.data.lapint{scalpN}, p.mesh.data.keepindex{scalpN}, p.mesh.data.repindex{scalpN}] = mesh_laplacian_interp(p.mesh.data.lapmat{scalpN}, elecindex);
end


fprintf('MESH_SCALP_INTERP: Multiplying Voltage by Interpolation matrix...\n');
tic;
% Interpolate voltage to all scalp mesh vertices (once off calc)
if isempty(p.mesh.data.repindex{scalpN}),
    p.mesh.data.Cdata{scalpN} = p.mesh.data.lapint{scalpN} * p.volt.data';
    %p.mesh.data.Cdata{scalpN} = p.volt.data * p.mesh.data.lapint{scalpN}';
else
    p.mesh.data.Cdata{scalpN} = p.mesh.data.lapint{scalpN} * p.volt.data(:,p.mesh.data.keepindex{scalpN})';
    %p.mesh.data.Cdata{scalpN} = p.volt.data(:,p.mesh.data.keepindex{scalpN}) * p.mesh.data.lapint{scalpN}';
end
t = toc;
fprintf('done (%6.2f sec).\n',t);


% Now a quick test of the interpolation: the values at the first
% Nelec vertices of the scalp mesh must be equal to the values
% of the original voltage data
TMP = p.mesh.data.Cdata{scalpN};
Nelec = size(p.volt.data,2);
if ~isequal(p.volt.data',TMP(1:Nelec,:)),
    msg = sprintf('MESH_SCALP_INTERP: Fatal Error in Mesh Voltage Interpolation, Reload all Data?\n');
    error(msg);
end


% At present, all scalp vertices are interpolated.  It might
% be useful at this point to determine all vertices that are
% "well below" the electrode array, say > 5cm from a nearest
% electrode and then set their voltage to zero.



return
