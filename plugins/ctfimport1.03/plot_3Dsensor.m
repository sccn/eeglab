function hf = plot_3Dsensor(channelfile,M,mesh,flat)
%function hf = plot_3Dsensor(channelfile,M,mesh,flat)
%
%Plots the channels, including the topography of measurements M, with or
%without a mesh
%
%inputs:
%channelfile: a file that contains a Channel structure, produced by BrainStorm
%M: a (nSensors x 1) vector with the measurements on the sensors. Use zeros
%    if you only want to see the sensors
%mesh: 'y' or 'n', plots a mesh around the sensors
%flat: 'y' or 'n', to get a flat mesh
%
%eg, use:
%plot_3Dsensor('GA_HiN_CueLInvNT-fEOG_channel.mat',zeros(275,1),'y')
%
% Author: Dimitrios Pantazis, PhD Student, USC, February 2005

%load channel file
load(channelfile)

%get channel locations
ndx = 1;
for i=1:size(Channel,2) %for all channels
    if strfind(Channel(i).Type,'MEG') & isempty(strfind(Channel(i).Type,'REF')) %if MEG
        MEGndx(ndx) = i;
        chanlocs(:,ndx) = Channel(i).Loc(:,1);
        ndx = ndx+1;
    end
end
chanlocs = chanlocs';    

%if M does not exist, assign zeros
if(~exist('M'))
    M  = zeros(size(MEGndx,2),1);
end

if(~exist('mesh','var'))
    mesh = 'y';
end

%create patch (Sylvain's code)
chanlocs(:,3)=chanlocs(:,3)-max(chanlocs(:,3));
[TH,PHI,R]=cart2sph(chanlocs(:,1),chanlocs(:,2),chanlocs(:,3));
%PHI2=zeros(size(PHI));
R2=R./cos(PHI).^.2;
[Y,X]=pol2cart(TH,R2);
%Y et X sont les coordonnées projetees dans le plan 'tangent superieur' pour chaque capteur
bord = convhull(Y,X);
ncapt = size(chanlocs,1);
[center,R] = bestfitsph(chanlocs);
coordC = chanlocs-(ones(ncapt,1)*center');
tri = convhulln(coordC./(norlig(coordC)'*ones(1,3)));
keep = find(~(ismember(tri(:,1),bord) & ismember(tri(:,2),bord) & ismember(tri(:,3),bord)));
tri = tri(keep,:);

%plot patch
hf = figure;
if exist('flat') & strcmp(flat,'y')
    h =trisurf(tri,X,Y,zeros(273,1),M); %flat
    view(90,270)
else
    h =trisurf(tri,chanlocs(:,1),chanlocs(:,2),chanlocs(:,3),M);
    view(0,90);
end
set(h,'edgecolor','none','facecolor','interp')
axis equal
axis off
axis vis3d



if(strcmp(mesh,'y'))
    set(h,'edgecolor','black');
    set(h,'marker','*');
end

set(gcf,'color','white');
