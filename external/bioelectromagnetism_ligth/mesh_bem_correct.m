function [p] = mesh_bem_correct(p,origin,dist,coin,cortex)

% mesh_bem_correct - Correct intersections of surfaces
% 
%[p] = mesh_bem_correct(p, origin, dist, coin, cortex)
% 
% p - the eeg_toolbox_defaults struct; in this case it
% is important that it contains p.mesh.data with four cells, 
% corresponding, in order, of cortex, inner skull,
% outer skull and scalp meshes.  If any of these surfaces
% are not available, the corresponding cell should be
% left empty, but allocated.  See mesh_bem_shells
% to obtain these surfaces.  The vertex coordinates 
% should be given in meters.
% 
% origin - 1x3, assumed to be [0 0 0] in Cartesian coordinates.
% 
% dist - the minimum distance, in mm, between each layer of
% the boundary element model (BEM).  The function converts
% these values from mm to meters.  At present, it assumes
% four layers in the BEM (modify the code otherwise).  For 
% example:
% 
%     dist.four2three = 6;    % scalp to outer skull
%     dist.four2two   = 8;    % scalp to inner skull
%     dist.three2two  = 5;    % outer skull to inner skull
%     dist.two2one    = 2;    % inner skull to cortex
% 
% coin  - this function assumes that each vertex in scalp and
% skull layers are coincident (which they are when created
% with mesh_bem_shells).  If these surface vertices are not
% coincident, set this input parameter to zero.  The inner 
% skull to cortex vertices are assumed not to be coincident.
% However, note, at present, this function will not work for
% surfaces that are not coincident (this option may be enabled
% in future development, hack the code otherwise).
% 
% cortex - 1 to process cortex/inner skull, 0 to process other
% surfaces only.  As the cortex/inner skull is the most time
% consuming process, this option allows for much faster 
% processing of the other surfaces.  This is useful for
% refined correction with adjustments of the dist struct.
% 
% There are no MRI intensity checks in this function, it
% simply adjusts the shells so there are no intersections.  The
% inner skull is first located outside the cortex (given cortex
% argument is 1). The inner skull is then static.  Next, the 
% scalp is located a given distance outside the inner skull
% (usually moved at the inferior parietal and occipito-temporal 
% regions, if at all).  The scalp is then static.  The outer
% skull is then located some distance inside the scalp and
% then located some distance outside the inner skull.  A final
% check and iterative refined movement ensures that the outer 
% skull is located between the scalp and the inner skull.  Note
% that the scalp to outer skull distance can be larger than
% the final placement of the outer skull, as the outer skull
% is finally relocated with respect to the inner skull.
% 

% $Revision: 1.1 $ $Date: 2009-04-28 22:13:56 $

% Licence:  GNU GPL, no implied or express warranties
% History:  10/2002, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nMESH_BEM_CORRECT: ...\n'); tic;


%fprintf('...sorry, undergoing repairs!\n\n'); return



if ~exist('origin','var'), origin = [0 0 0];
elseif isempty(origin),    origin = [0 0 0];
end
xo = origin(1); yo = origin(2); zo = origin(3);

if ~exist('dist','var'),
    dist.four2three = 6;    % scalp to outer skull
    dist.four2two   = 8;    % scalp to inner skull
    dist.three2two  = 5;    % outer skull to inner skull
    dist.two2one    = 2;    % inner skull to cortex
elseif isempty(dist),
    dist.four2three = 6;    % scalp to outer skull
    dist.four2two   = 8;    % scalp to inner skull
    dist.three2two  = 5;    % outer skull to inner skull
    dist.two2one    = 2;    % inner skull to cortex
end

% convert dist from mm to meters
dist.four2three = dist.four2three / 1000;
dist.four2two   = dist.four2two   / 1000;
dist.three2two  = dist.three2two  / 1000;
dist.two2one    = dist.two2one    / 1000;


if ~exist('coin','var'),
    coin = 1;
elseif isempty(dist),
    coin = 1;
end

if ~exist('cortex','var'),
    % run correction for inner skull/cortex
    cortex = 1;
elseif isempty(cortex),
    cortex = 1;
end

% Identify each layer of the BEM
[scalpN,exists]  = mesh_check(p,'scalp');
if ~exists, error('...No scalp in p struct');
else,       fprintf('...found cell %d of p.mesh.data contains scalp\n',scalpN);
end
[oskullN,exists] = mesh_check(p,'outer_skull');
if ~exists, error('...No outer_skull in p struct');
else,       fprintf('...found cell %d of p.mesh.data contains outer_skull\n',oskullN);
end
[iskullN,exists] = mesh_check(p,'inner_skull');
if ~exists, error('...No inner_skull in p struct');
else,       fprintf('...found cell %d of p.mesh.data contains inner_skull\n',iskullN);
end
[cortexN,exists] = mesh_check(p,'cortex');
if ~exists,
    [cortexN,exists] = mesh_check(p,'pial');
    if ~exists, error('...No cortex in p struct');
    else fprintf('...found cell %d of p.mesh.data contains cortex\n',cortexN);
    end
else,    fprintf('...found cell %d of p.mesh.data contains cortex\n',cortexN);
end


% Find distances and direction cosines of all surface vertices
fprintf('...calculating vertex radius and direction cosines: scalp\n');
Xscalp = p.mesh.data.vertices{scalpN}(:,1) - xo;
Yscalp = p.mesh.data.vertices{scalpN}(:,2) - yo;
Zscalp = p.mesh.data.vertices{scalpN}(:,3) - zo;
Dscalp = sqrt( (Xscalp).^2 + (Yscalp).^2 + (Zscalp).^2 );
Lscalp = (Xscalp)./Dscalp; % cos alpha
Mscalp = (Yscalp)./Dscalp; % cos beta
Nscalp = (Zscalp)./Dscalp; % cos gamma

fprintf('...calculating vertex radius and direction cosines: outer skull\n');
Xoskull = p.mesh.data.vertices{oskullN}(:,1) - xo;
Yoskull = p.mesh.data.vertices{oskullN}(:,2) - yo;
Zoskull = p.mesh.data.vertices{oskullN}(:,3) - zo;
Doskull = sqrt( (Xoskull).^2 + (Yoskull).^2 + (Zoskull).^2 );
Loskull = (Xoskull)./Doskull; % cos alpha
Moskull = (Yoskull)./Doskull; % cos beta
Noskull = (Zoskull)./Doskull; % cos gamma

fprintf('...calculating vertex radius and direction cosines: inner skull\n');
Xiskull = p.mesh.data.vertices{iskullN}(:,1) - xo;
Yiskull = p.mesh.data.vertices{iskullN}(:,2) - yo;
Ziskull = p.mesh.data.vertices{iskullN}(:,3) - zo;
Diskull = sqrt( (Xiskull).^2 + (Yiskull).^2 + (Ziskull).^2 );
Liskull = (Xiskull)./Diskull; % cos alpha
Miskull = (Yiskull)./Diskull; % cos beta
Niskull = (Ziskull)./Diskull; % cos gamma

% Note that the L,M,N direction cosines of the above
% should be equal when these surface vertices are
% coincident (given 1-2 decimal places).

% check that max distances are less than 1 meter
if max(Dscalp) > 1,  error('...maximum scalp radius > 1 meter.\n'); end
if max(Doskull) > 1, error('...maximum outer skull radius > 1 meter.\n'); end
if max(Diskull) > 1, error('...maximum inner skull radius > 1 meter.\n'); end


CORRECTION_ADJUSTMENT = 0.0001;


if cortex,
    
    iSkull.vertices = p.mesh.data.vertices{iskullN};
    iSkull.faces    = p.mesh.data.faces{iskullN};
    
    cortex.vertices = p.mesh.data.vertices{cortexN};
    cortex.faces    = p.mesh.data.faces{cortexN};
    
    fprintf('...processing cortex/inner skull (a lengthy process)...\n');
    fprintf('...tesselating cortex vertices for searching...'); tic;
    tri = delaunayn(cortex.vertices);
    t = toc; fprintf(' (%5.2f sec)\n',t);
    
    correct = 1;
    while correct,
        % Find cortex vertices closest to the inner skull vertex
        iSkullNvert = size(iSkull.vertices,1);
        fprintf('...checking for nearest cortex vertex to %d inner skull vertices...',iSkullNvert); tic;
        cortex_vertex_index = dsearchn(cortex.vertices,tri,iSkull.vertices);
        t = toc; fprintf(' (%5.2f sec)\n',t);
        
        Xcortex = cortex.vertices(cortex_vertex_index,1) - xo;
        Ycortex = cortex.vertices(cortex_vertex_index,2) - yo;
        Zcortex = cortex.vertices(cortex_vertex_index,3) - zo;
        Dcortex = sqrt( (Xcortex).^2 + (Ycortex).^2 + (Zcortex).^2 );
        if max(Dcortex) > 1, error('...maximum cortex radius > 1 meter.\n'); end
        
        % Don't calculate direction cosines for cortex, because these
        % vertices are not to be moved
        
        % ensure the inner skull vertices are outside the cortex
        dif = Diskull - Dcortex;
        correct = find(dif < dist.two2one);
        
        if correct,
            
            % plot the cortex and inner skull
            iSkullfig = figure;
            patch('vertices',cortex.vertices,'faces',cortex.faces,'facecolor',[.7 0  0],'edgecolor','none'); hold on
            patch('vertices',iSkull.vertices,'faces',iSkull.faces,'facecolor',[ 0 0 .7],'edgecolor','none','facealpha',.4);
            axis vis3d; view(3); axis tight; light
            drawnow
            
            vert = iSkull.vertices(correct,:);
            x = vert(:,1);
            y = vert(:,2);
            z = vert(:,3);
            plot3(x,y,z,'bo')
            
            pause(5)
            close(iSkullfig);
            
            fprintf('...correcting %4d inner skull vertices < cortex + %d mm\n',size(correct,1),dist.two2one * 1000);
            Diskull(correct) = Diskull(correct) + (dist.two2one - dif(correct)) + CORRECTION_ADJUSTMENT;
            Xiskull = (Liskull .* Diskull);
            Yiskull = (Miskull .* Diskull);
            Ziskull = (Niskull .* Diskull);
            iSkull.vertices = [ Xiskull Yiskull Ziskull ];
            iSkull = mesh_smooth(iSkull,origin);
            
        else
            fprintf('...all inner skull vertices > cortex + %d mm\n',dist.two2one * 1000);
        end
        % iterate until corrected, finding cortex vertices
        % each time, as the nearest cortex vertex can change
        if ishandle(iSkullfig), close(iSkullfig); end
    end
    
    p.mesh.data.vertices{iskullN} = iSkull.vertices;
    p.mesh.data.faces{iskullN}    = iSkull.faces;
    clear iSkull cortex;
    
end

% ensure the scalp vertices are outside the inner skull
FV.vertices = p.mesh.data.vertices{scalpN};
FV.faces    = p.mesh.data.faces{scalpN};

correct = 1;
while correct,
    dif = Dscalp - Diskull;
    correct = find(dif < dist.four2two);
    if correct,
        fprintf('...correcting %4d scalp vertices < inner skull + %d mm\n',size(correct,1),dist.four2two * 1000);
        Dscalp(correct) = Dscalp(correct) + (dist.four2two - dif(correct)) + CORRECTION_ADJUSTMENT;
        Xscalp = (Lscalp .* Dscalp);
        Yscalp = (Mscalp .* Dscalp);
        Zscalp = (Nscalp .* Dscalp);
        FV.vertices = [ Xscalp Yscalp Zscalp ];
        FV = mesh_smooth(FV,origin);
    else
        fprintf('...all scalp vertices > inner skull + %d mm\n',dist.four2two * 1000);
    end
end
p.mesh.data.vertices{scalpN} = FV.vertices;
p.mesh.data.faces{scalpN}    = FV.faces;
clear FV;

% ensure the outer skull vertices are inside the scalp
FV.vertices = p.mesh.data.vertices{oskullN};
FV.faces    = p.mesh.data.faces{oskullN};

correct = 1;
while correct,
    dif = Dscalp - Doskull;
    correct = find(dif < dist.four2three);
    if correct,
        fprintf('...correcting %4d outer skull vertices > scalp - %d mm\n',size(correct,1),dist.four2three * 1000);
        Doskull(correct) = Doskull(correct) - (dist.four2three - dif(correct)) - CORRECTION_ADJUSTMENT;
        Xoskull = (Loskull .* Doskull);
        Yoskull = (Moskull .* Doskull);
        Zoskull = (Noskull .* Doskull);
        FV.vertices = [ Xoskull Yoskull Zoskull ];
        FV = mesh_smooth(FV,origin);
    else
        fprintf('...all outer skull vertices < scalp - %d mm\n',dist.four2three);
    end
end

% ensure the outer skull vertices are outside the inner skull
correct = 1;
while correct,
    dif = Doskull - Diskull;
    correct = find(dif < dist.three2two);
    if correct,
        fprintf('...correcting %4d outer skull vertices < inner skull + %d mm\n',size(correct,1),dist.three2two * 1000);
        Doskull(correct) = Doskull(correct) + (dist.three2two - dif(correct)) + CORRECTION_ADJUSTMENT;
        Xoskull = (Loskull .* Doskull);
        Yoskull = (Moskull .* Doskull);
        Zoskull = (Noskull .* Doskull);
        FV.vertices = [ Xoskull Yoskull Zoskull ];
        FV = mesh_smooth(FV,origin);
    else
        fprintf('...all outer skull vertices > inner skull + %d mm\n',dist.three2two * 1000);
    end
end


% finally, ensure the outer skull vertices are between the scalp and the inner skull
measure = .002; %start with 2 mm
pass = 0;
move = [1 1];
while find(move),
    
    pass = pass + 1;
    if pass > 1,
        % decrease the testing measure, to allow progressive
        % determination of outer skull within scalp/inner skull
        measure = measure - 0.00025;
    end
    
    % Confirm that oskull is at least 2mm inside scalp
    scalp2oskull = Dscalp - Doskull;
    correct = find(scalp2oskull < measure);
    if correct,
        fprintf('...correcting %4d outer skull vertices > scalp - %d mm\n',size(correct,1),measure * 1000);
        Doskull(correct) = Doskull(correct) - (measure - scalp2oskull(correct)) - CORRECTION_ADJUSTMENT;
        Xoskull = (Loskull .* Doskull);
        Yoskull = (Moskull .* Doskull);
        Zoskull = (Noskull .* Doskull);
        FV.vertices = [ Xoskull Yoskull Zoskull ];
        move(1) = 1;
    else,
        fprintf('...confirmed outer skull vertices < scalp - %d mm\n',measure * 1000);
        move(1) = 0;
    end
    
    oskull2iskull = Doskull - Diskull;
    correct = find(oskull2iskull < measure);
    if correct,
        fprintf('...correcting %4d outer skull vertices < inner skull + %d mm\n',size(correct,1),measure * 1000);
        Doskull(correct) = Doskull(correct) + (measure - oskull2iskull(correct)) + CORRECTION_ADJUSTMENT;
        Xoskull = (Loskull .* Doskull);
        Yoskull = (Moskull .* Doskull);
        Zoskull = (Noskull .* Doskull);
        FV.vertices = [ Xoskull Yoskull Zoskull ];
        move(2) = 1;
    else,
        fprintf('...confirmed outer skull vertices > inner skull + %d mm\n',measure * 1000);
        move(2) = 0;
    end
end
p.mesh.data.vertices{oskullN} = FV.vertices;
p.mesh.data.faces{oskullN}    = FV.faces;
clear FV;


t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code can be used to plot the vertices that are to be
% corrected
% 
%     v = correct;
%     index = [scalpN oskullN];
%     
%     patch('vertices',p.mesh.data.vertices{index(1)},'faces',p.mesh.data.faces{index(1)},...
%           'FaceColor',[.9 .9 .9],'Edgecolor','none','FaceAlpha',.6);
%     daspect([1 1 1]); axis tight; hold on
%     
%     vert = p.mesh.data.vertices{index(1)};
%     x = vert(v,1);
%     y = vert(v,2);
%     z = vert(v,3);
%     plot3(x,y,z,'ro')
%     
%     patch('vertices',p.mesh.data.vertices{index(2)},'faces',p.mesh.data.faces{index(2)},...
%           'FaceColor',[.0 .6 .0],'Edgecolor',[.2 .2 .2],'FaceAlpha',.8);
%     
%     vert = p.mesh.data.vertices{index(2)};
%     x = vert(v,1);
%     y = vert(v,2);
%     z = vert(v,3);
%     plot3(x,y,z,'bo')
