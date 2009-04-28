function[p] = elec_load_brainstorm(p)

% elec_load_brainstorm - Load brainstorm Channel data into eeg_toolbox
% 
% Usage: [p] = elec_load_brainstorm(p)
% 
% where:    p.elec.path & p.elec.file define the brainsform file
% 
% The brainstorm file consists of Channel, which is an array of 
% structures.  The fields are:
%           
%     Loc     - a 3x2 matrix of electrode coordinates (x,y,z in rows).
%               BrainStorm (x,y,z meters) = 3Dspace (x,y,z cm) / 100.
%     Orient  - a corresponding matrix of sensor orientations (MEG); 
%               all zero for EEG.
%     Weight  - a vector of relative or absolute weights (eg, amplification);
%               all ones for this routine.
%     Type    - a character string, 'EEG' in this case.
%     Name    - a charater string indicating the electrode name.
%     Comment - a charater string indicating the reference electrode.  Empty
%               for active electrodes and 'EEG REF' for reference electrode.
% 
% See brainstorm website at http://neuroimage.usc.edu/, including a
% download pdf file describing the brainstorm database formats.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no express or implied warranties
% History:  11/2002, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('p', 'var'), error('...no p struct input\n');
elseif isempty(p),     error('...no p struct input\n');
end

fprintf('\nELEC_LOAD_BRAINSTORM...\n'); tic;

fprintf('...Converting brainstorm to p structure.\n');


[path,file,ext] = fileparts(strcat(p.elec.path,filesep,p.elec.file));
bsfile = fullfile(path,[file ext]);

load(bsfile, 'Channel');

for i=1:length(Channel),
    
    if findstr(lower(Channel(i).Comment),'ref'),
        p.elec.data.ref = Channel(i).Loc(:,1)';
    else
        p.elec.data.x(i,1) = Channel(i).Loc(1,1);
        p.elec.data.y(i,1) = Channel(i).Loc(2,1);
        p.elec.data.z(i,1) = Channel(i).Loc(3,1);
        p.elec.data.label{i,1} = Channel(i).Name;
    end
end

p.elec.data.centroid = [0 0 0];
p.elec.data.nasion = [];
p.elec.data.lpa = [];
p.elec.data.rpa = [];

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return
