
% eeg_colormap_script - play with colour mapping

clear all, close all

Map = eeg_colormap('Red/Blue/White'); colormap(Map);

for i = 1:size(Map,1)

    c = Map(i,:);   % define patch color

    xp = [-5  5  5 -5];
    yp = [ i  i  i+1 i+1];
    
    patch(xp,yp,c)

end
colorbar
