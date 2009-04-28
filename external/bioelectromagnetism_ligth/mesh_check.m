function [n,exists] = mesh_check(p,meshtype)

% mesh_check - Check for mesh data in p struct
%
% Usage: [n,exists] = mesh_check(p,meshtype)
%
% p is the eeg_toolbox struct (see eeg_toolbox_defaults)
% and meshtype is a string argument that is compared
% with all cells of p.mesh.data.meshtype.  If it matches,
% the return value 'n' is the matching cell, ie,
% p.mesh.data.meshtype{n} and the boolean 'exists'
% variable is true. This function is called from
% mesh_open and other mesh utilities to facilitate
% efficient handling of the p.mesh.data struct.
%

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:57 $

% Licence:  GNU GPL, no implied or express warranties
% History:  02/2002, Darren.Weber_at_radiology.ucsf.edu
%           02/2004, Darren.Weber_at_radiology.ucsf.edu
%                    check if the current cell is empty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 1;          % this is the first mesh in p.mesh.data
exists = 0;

if ~isempty(p.mesh.data),
    if isfield(p.mesh.data,'meshtype'),
        
        Nmeshes = length(p.mesh.data.meshtype);
        
        % new mesh can be added at the end
        n = Nmeshes + 1;
        
        % check if the p.mesh.current is empty
        if p.mesh.current <= Nmeshes,
          if isempty(p.mesh.data.meshtype{p.mesh.current}),
            n = p.mesh.current;
          end
        end
        
        % check whether any other cells contain the same mesh data so it
        % can be replaced
        indices = strcmp(p.mesh.data.meshtype,meshtype);
        found = find(indices);
        if found,
          n = found(1);
          exists = 1;
        end
        
    end
end

return
