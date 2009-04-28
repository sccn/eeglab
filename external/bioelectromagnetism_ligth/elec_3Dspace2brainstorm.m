function [Channel] = elec_3dspace2brainstorm(filename,N,bsfile)

% elec_3dspace2brainstorm - Convert NeuroScan 3Dspace ascii to brainstorm file
% 
% The ascii file format is that of NeuroScan 3Dspace export files.  Each row
% of the file comprises an electrode label, an electrode type code, and
% the x,y,z coordinates (cm).  Each field is separated by spaces.  See
% ELEC_LOAD for more information.
% 
% Useage: Channel = elec_3dspace2brainstorm(filename,[N],[brainstormfile])
% 
% where:    file = 'path\filename' with format described in ELEC_LOAD
% 
%           N is how many electrodes to load (rows of 3Dspace file, 129 default).
% 
% Result:   Channel is an array of structures.  The fields are:
%           
%           Loc     - a 3x2 matrix of electrode coordinates (x,y,z in rows).
%                     BrainStorm (x,y,z meters) = 3Dspace (x,y,z cm) / 100.
%           Orient  - a corresponding matrix of sensor orientations (MEG); 
%                     all zero for EEG.
%           Weight  - a vector of relative or absolute weights (eg, amplification);
%                     all ones for this routine.
%           Type    - a character string, 'EEG' in this case.
%           Name    - a charater string indicating the electrode name.
%           Comment - a charater string indicating the reference electrode.  Empty
%                     for active electrodes and 'EEG REF' for reference electrode.
% 
% See brainstorm website at http://neuroimage.usc.edu/, including a
% download pdf file describing the brainstorm database formats.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:54 $

% Licence:  GNU GPL, no express or implied warranties
% History:  20/05/2002, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('N', 'var'), N = 129;
elseif isempty(N),     N = 129;
end

file = char(filename);

[elec,type,X,Y,Z] = elec_load(file,[],[],[],[],N);

tic;

fprintf('\nELEC_3DSPACE2BRAINSTORM...\n');
fprintf('...Converting to brainstorm structure.\n');

elecindex = find(type == 69);
Ee = elec(elecindex);
Xe = X(elecindex) ./ 100;   % 3Dspace is cm, BrainStorm is m
Ye = Y(elecindex) ./ 100;
Ze = Z(elecindex) ./ 100;

reftype = ones(size(type)) * 120;
refindex = find(type == reftype);
Eref = elec(refindex);
Xref = X(refindex) ./ 100;   % 3Dspace is cm, BrainStorm is m
Yref = Y(refindex) ./ 100;
Zref = Z(refindex) ./ 100;

for i=1:length(elecindex),
    Channel(i).Loc = [[Xe(i) Ye(i) Ze(i)]',[Xref Yref Zref]'];
    Channel(i).Orient = [];     % used for MEG rather than EEG
    Channel(i).Weight = 1;      % Like Amplification
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

if ~exist('bsfile', 'var'),
    bsfile = 'channel';
elseif isempty(bsfile),
    bsfile = 'channel';
end

if findstr('.mat',bsfile),
    bsfile = strrep(bsfile,'.mat','');
end

fprintf('...saving BrainStorm channel data to:\n\t%s.mat\n',bsfile);
save(bsfile, 'Channel');


t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return
