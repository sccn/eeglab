function [p] = elec_write_brainstorm(p)

% elec_write_brainstorm - Write p.elec.data to brainstorm file
% 
% Useage: [p] = elec_write_brainstorm(p)
% 
% The electrode orientation saved is consistent with the BrainStorm
% 'Neuromag' orientation system, where +x is through the right ear
% and +y is through the nasion.
% 
% See brainstorm website at http://neuroimage.usc.edu/, including a
% download pdf file describing the brainstorm database formats.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:55 $

% Licence:  GNU GPL, no express or implied warranties
% History:  10/2002, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

fprintf('\nELEC_WRITE_BRAINSTORM...\n');
fprintf('...Converting to brainstorm structure.\n');

Ee = p.elec.data.label;
Xe = p.elec.data.x;
Ye = p.elec.data.y;
Ze = p.elec.data.z;

Eref = 'Ref';
Xref = p.elec.data.ref(1);
Yref = p.elec.data.ref(2);
Zref = p.elec.data.ref(3);

for i=1:length(Xe),
    Channel(i).Loc = [[Xe(i) Ye(i) Ze(i)]',[Xref Yref Zref]']; % See note below
    Channel(i).Orient = [];     % used for MEG rather than EEG
    Channel(i).Weight = 1;      % used for MEG rather than EEG
    Channel(i).Type = 'EEG';
    Channel(i).Name = char(Ee(i));
    Channel(i).Comment = '';
end
Channel(i+1).Loc = [[Xref Yref Zref]',[Xref Yref Zref]'];
Channel(i+1).Orient = [];
Channel(i+1).Weight = 1;
Channel(i+1).Type = 'EEG';
Channel(i+1).Name = char(Eref);
Channel(i+1).Comment = 'EEG REF';


[path,name,ext] = fileparts(strcat(p.elec.path,filesep,p.elec.file));
name = sprintf('%s_channel',name);
ext  = '.mat';
file = fullfile(path,[name ext]);

fprintf('...saving BrainStorm channel data to:\n\t%s\n',file);
save(file, 'Channel');

% channel(i).Loc field replicates the X,Y,Z coordinates of each
% electrode to conform with a format that accommodates MEG coordinates,
% which include two sensors per channel.

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return
