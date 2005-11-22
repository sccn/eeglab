% coregister() -  Coregister electrode location file with template or mesh.
%
% Usage:
%        >> coregister( chan1, chan2, 'key', 'val' )
%
% Inputs:
%    chan1    - channel location file or eeglab channel location 
%               structure to align.
%    chan2    - reference channel location file or eeglab channel location 
%               structure
%
% Optional input:
%    'alignfid'  - [cell array of string] name of fiducials for alignment.
%    'warp'      - [cell array of string] name of electrode for warping.
%    'transform' - [real array] homogenous transformation matrix or 1x9
%                  matrix containing 9 parameter traditional "Talairach-model" 
%                  transformation (see traditional()).
%    'mesh'      - [cell array|string] head mesh. Can contain 
%                  { points triangles } or { point triangle normals } see the 
%                  function plotmesh() for details. Can also contain the name
%                  of a file containing head mesh information (several format
%                  recognized).
%    'autoscale' - ['on'|'off'] autoscale electrode radius when aligning 
%                  fiducials default is 'on'.
%
% Output:
%    chan1       - transformed channel location structure
%    transform   - transformation matrix. Use function traditional() to 
%                  convert to homogenous transformation matrix and input
%                  it into functions like headplot().
% 
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2005
        
%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, SCCN, INC, UCSD, 2005, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: not supported by cvs2svn $
% Revision 1.2  2005/11/08 23:08:55  arno
% header
%
% Revision 1.1  2005/11/08 22:48:05  arno
% Initial revision
%

function [ chan1, transformmat ] = coregister(chan1, chan2, varargin)

if nargin < 1
    % chan1 = readlocs('/home/arno/P3FigArno/juliedata/jop3_raw.elp');
    chan1 = readlocs('D:\data\arnodata\ad-256.elp');
   % manually for Arno:  transf = [4 0 -50 -0.3 0 -1.53 1.05 1.1 1.1]
   % manually for Julie: transf = [-4 -6 -50 -0.37 0 -1.35 1.1 1.15 1.1]
    chan2 = readlocs('D:\matlab\eeglab\plugins\dipfit2.0\standard_BEM\elec\standard_1005.elc');
    normalize = 1;
end;

% internal command
% ----------------
if isstr(chan1) 
    com = chan1;
    fid = chan2;

    % update GUI
    % ----------
    if strcmpi(com, 'redraw'), redrawgui(fid); return; end;

    % select electrodes and warp montage
    % ----------------------------------
    dat = get(fid, 'userdata');
    if strcmpi(com, 'fiducials')
        [clist1 clist2] = pop_chancoresp( dat.elec1, dat.elec2, 'autoselect', 'fiducials', 'gui', 'off');
        [ tmp dat.transform ] = align_fiducials(dat.elec1, dat.elec2, dat.elec2.label(clist2));
    elseif strcmpi(com, 'warp')
        [clist1 clist2] = pop_chancoresp( dat.elec1, dat.elec2, 'autoselect', 'fiducials');
        % copy electrode names
        tmpelec2 = dat.elec2;
        for index = 1:length(clist2)
            tmpelec2.label{clist2(index)} = dat.elec1.label{clist1(index)};
        end;
        [ tmp dat.transform ] = warp_chans(dat.elec1, dat.elec2, tmpelec2.label(clist2));
    end;
    set(fid, 'userdata', dat);
    redrawgui(fid); 
    return;
    
end;

% check input arguments
% ---------------------
defaultmesh = 'D:\matlab\eeglab\plugins\dipfit2.0\standard_BEM\standard_vol.mat';
g = finputcheck(varargin, { 'alignfid'   'cell'  {}      {};
                            'warp'       'cell'  {}      {};
                            'transform'  'real'  []      [];
                            'autoscale'  'string' { 'on' 'off' } 'on';
                            'mesh'       ''      []   defaultmesh });
if isstr(g), error(g); end;

% load mesh if any
% ----------------
if ~isempty(g.mesh)
    if isstr(g.mesh)
        g.mesh  = load(g.mesh);
    end;
    if isstruct(g.mesh)
        if isfield(g.mesh, 'vol')
            dat.meshpnt = g.mesh.vol.bnd(1).pnt;
            dat.meshtri = g.mesh.vol.bnd(1).tri;
        elseif isfield(g.mesh, 'bnd')
            dat.meshpnt = g.mesh.bnd(1).pnt;
            dat.meshtri = g.mesh.bnd(1).tri;
        elseif isfield(g.mesh, 'TRI1')
            dat.meshpnt = g.mesh.POS;
            dat.meshtri = g.mesh.TRI1;
        else
            error('Unknown Matlab mesh file');
        end;
    else
        dat.meshpnt = g.mesh{1};
        dat.meshtri = g.mesh{2};
    end;
else
    dat.meshpnt = [];
    dat.meshtri = [];
end;

% transform to arrays chan1
% -------------------------
[tmp1 tmp2 tmp3 ind1] = readlocs(chan1);
TMP   = eeg_emptyset;
TMP.chanlocs = chan1;
cfg   = eeglab2fieldtrip(TMP, 'timelockanalysis');
elec1 = cfg.elec;

% transform to arrays chan2
% -------------------------
if ~isempty(chan2)
    TMP.chanlocs = chan2;
    cfg = eeglab2fieldtrip(TMP, 'timelockanalysis');
    elec2 = cfg.elec;
else 
    elec2 = [];
    dat.transform = [ 0 0 0 0 0 0 1 1 1 ];
end;

% copy or compute alignment matrix
% --------------------------------
if ~isempty(g.transform)
    dat.transform = g.transform;
else

    % perfrom alignment
    % -----------------
    if ~isempty(g.alignfid)
        % autoscale
        % ---------
        electmp = elec1;
        if strcmpi(g.autoscale, 'on')
            avgrad1 = sqrt(sum(elec1.pnt.^2,2));
            avgrad2 = sqrt(sum(elec2.pnt.^2,2));
            ratio   = mean(avgrad2)/mean(avgrad1)*1.05;
            transmat = [ratio 0 0 0; 0 ratio 0 0; 0 0 ratio 0; 0 0 0 1];
            electmp.pnt = elec1.pnt*ratio;
        end;

        [ electransf dat.transform ] = align_fiducials(electmp, elec2, g.alignfid);
        %dat.transform = invtrad(dat.transform, 1)
        
        if strcmpi(g.autoscale, 'on')
            dat.transform =dat.transform*transmat;
        end;        
        
    elseif ~isempty(g.warp)
        [ electransf dat.transform ] = warp_chans(elec1, elec2, g.warp);
    else
        dat.transform = [0 0 0 0 0 0 1 1 1];
    end;
    
end;

% find common electrode names
% ---------------------------
dat.elec1      = elec1;
dat.elec2      = elec2;
fid = figure('userdata', dat);

if 1
    header    = 'dattmp = get(gcbf, ''userdata'');';
    footer    = 'set(gcbf, ''userdata'', dattmp); clear dattmp; coregister(''redraw'', gcbf);';
    cbpitch   = [ header 'dattmp.transform(4) = str2num(get(gcbo, ''string''));' footer ];
    cbroll    = [ header 'dattmp.transform(5) = str2num(get(gcbo, ''string''));' footer ];
    cbyaw     = [ header 'dattmp.transform(6) = str2num(get(gcbo, ''string''));' footer ];
    cbresizex = [ header 'dattmp.transform(7) = str2num(get(gcbo, ''string''));' footer ];
    cbresizey = [ header 'dattmp.transform(8) = str2num(get(gcbo, ''string''));' footer ];
    cbresizez = [ header 'dattmp.transform(9) = str2num(get(gcbo, ''string''));' footer ];
    cbright   = [ header 'dattmp.transform(1) = str2num(get(gcbo, ''string''));' footer ];
    cbforward = [ header 'dattmp.transform(2) = str2num(get(gcbo, ''string''));' footer ];
    cbup      = [ header 'dattmp.transform(3) = str2num(get(gcbo, ''string''));' footer ];
    cb_ok     = 'set(gcbo, ''userdata'', ''ok'')';
    cb_warp   = 'coregister(''warp'', gcbf);';
    cb_fid    = 'coregister(''fiducials'', gcbf);';
    
    opt = { 'unit', 'normalized', 'position' };
    h = uicontrol( opt{:}, [0    .15  1  .02], 'style', 'text', 'string', '');
    h = uicontrol( opt{:}, [0    .1  .15 .05], 'style', 'text', 'string', 'Pitch (rad)');
    h = uicontrol( opt{:}, [0    .05 .15 .05], 'style', 'text', 'string', 'Roll (rad)' );
    h = uicontrol( opt{:}, [0    0   .15 .05], 'style', 'text', 'string', 'Yaw (rad)');
    h = uicontrol( opt{:}, [0.15 .1  .1  .05], 'tag', 'pitch', 'callback', cbpitch, 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.15 .05 .1  .05], 'tag', 'roll' , 'callback', cbroll , 'style', 'edit', 'string', '' );
    h = uicontrol( opt{:}, [0.15 0   .1  .05], 'tag', 'yaw'  , 'callback', cbyaw  , 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.25 .1  .15 .05], 'style', 'text', 'string', 'Resize {x}');
    h = uicontrol( opt{:}, [0.25 .05 .15 .05], 'style', 'text', 'string', 'Resize {y}' );
    h = uicontrol( opt{:}, [0.25 0   .15 .05], 'style', 'text', 'string', 'Resize {z}');
    h = uicontrol( opt{:}, [0.4  .1  .1  .05], 'tag', 'resizex', 'callback', cbresizex, 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.4  .05 .1  .05], 'tag', 'resizey', 'callback', cbresizey, 'style', 'edit', 'string', '' );
    h = uicontrol( opt{:}, [0.4  0   .1  .05], 'tag', 'resizez', 'callback', cbresizez, 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.5  .1  .2  .05], 'style', 'text', 'string', 'Right {x in mm}');
    h = uicontrol( opt{:}, [0.5  .05 .2  .05], 'style', 'text', 'string', 'Forward {y in mm}' );
    h = uicontrol( opt{:}, [0.5  0   .2  .05], 'style', 'text', 'string', 'Up {z in mm}');
    h = uicontrol( opt{:}, [0.7  .1  .1  .05], 'tag', 'right'  , 'callback', cbright  , 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.7  .05 .1  .05], 'tag', 'forward', 'callback', cbforward, 'style', 'edit', 'string', '' );
    h = uicontrol( opt{:}, [0.7  0   .1  .05], 'tag', 'up'     , 'callback', cbup     , 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.8  .1  .2  .05], 'style', 'pushbutton', 'string', 'Align fiducials', 'callback', cb_fid);
    h = uicontrol( opt{:}, [0.8  .05 .2  .05], 'style', 'pushbutton', 'string', 'Warp montage', 'callback', cb_warp );
    h = uicontrol( opt{:}, [0.8  0   .1  .05], 'style', 'pushbutton', 'string', 'Cancel', 'callback', 'close(gcbf);' );
    h = uicontrol( opt{:}, [0.9  0   .1  .05], 'style', 'pushbutton', 'string', 'Ok', 'tag', 'ok', 'callback', cb_ok);
end;

coregister('redraw', fid);
try, icadefs; set(gcf, 'color', BACKCOLOR); catch, end

% wait until button press and return values
% -----------------------------------------
waitfor( findobj('parent', fid, 'tag', 'ok'), 'userdata');
try, findobj(fid); % figure still exist ?
catch, transformmat = []; chan1 = []; return; end;
dat = get(fid, 'userdata');
transformmat = dat.transform;
chan1        = dat.electransf;
close(fid);

% plot electrodes
% ---------------
function plotelec(elec, color);
    
    X1 = elec.pnt(:,1);
    Y1 = elec.pnt(:,2);
    Z1 = elec.pnt(:,3);    

    XL = xlim;
    YL = ylim;
    ZL = zlim;
    
    lim=max(1.05*max([X1;Y1;Z1]), max([XL YL ZL]));
    eps=lim/20;
    delete(findobj(gcf, 'tag', [ color 'plot']));
    h = plot3(X1,Y1,Z1, [ color 'o' ]); hold on;
    set(h, 'tag', [ color 'plot'], ...
           'marker', '.', 'markersize', 20);

    
    %if ~isempty(elocname)
    %    plot3(X1(elec_ind),Y1(elec_ind),Z1(elec_ind),'b*')
    %end;

    % plot axis and labels
    %- -------------------
    plot3([0.08 0.12],[0 0],[0 0],'r','LineWidth',4) % nose
    plot3([0 lim],[0 0],[0 0],'b--')                 % axis
    plot3([0 0],[0 lim],[0 0],'g--')
    plot3([0 0],[0 0],[0 lim],'r--')
    plot3(0,0,0,'b+')
    text(lim+eps,0,0,'X','HorizontalAlignment','center',...
         'VerticalAlignment','middle','Color',[0 0 0],...
         'FontSize',10)
    text(0,lim+eps,0,'Y','HorizontalAlignment','center',...
         'VerticalAlignment','middle','Color',[0 0 0],...
         'FontSize',10)
    text(0,0,lim+eps,'Z','HorizontalAlignment','center',...
         'VerticalAlignment','middle','Color',[0 0 0],...
         'FontSize',10)
    box on
    
    %if ~isempty(elocname)
    %    for i = 1:length(elec_ind)
    %        text(X(elec_ind(i)),Y(elec_ind(i)),Z(elec_ind(i))+eps,elocname(elec_ind(i)),'HorizontalAlignment','center',...
    %             'VerticalAlignment','middle','Color',[0 0 0], 'FontSize',10)
    %    end
    %nd;

    axis([-lim lim -lim lim -lim*0.5 lim])
    axis equal;

% align fiducials
% ---------------
function [elec1, transf] = align_fiducials(elec1, elec2, fidnames)

    % rename fiducials in template
    % ----------------------------
    fidtemplate = { 'Nz' 'LPA' 'RPA' };
    ind = strmatch(fidtemplate{1}, elec2.label, 'exact'); elec2.label{ind} = fidnames{1};
    ind = strmatch(fidtemplate{2}, elec2.label, 'exact'); elec2.label{ind} = fidnames{2};
    ind = strmatch(fidtemplate{3}, elec2.label, 'exact'); elec2.label{ind} = fidnames{3};

    cfg          = [];
    cfg.elec     = elec1;
    cfg.template = elec2;
    cfg.method   = 'realignfiducial'; 
    cfg.fiducial = fidnames;
    elec3 = electrodenormalize(cfg);
    transf = elec3.m
    %transf(1,1:3) = transf(1,1:3)*1.23;
    %transf(2,1:3) = transf(2,1:3)*1.05;
    %transf(3,1:3) = transf(3,1:3)*1.05;
    elec1.pnt = transf*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
    transfmat = traditional(elec3.m);
    elec1.pnt = elec1.pnt(1:3,:)';

% warp channels
% -------------
function [elec1, transf] = warp_chans(elec1, elec2, chanlist)
    cfg          = [];
    cfg.elec     = elec1;
    cfg.template = elec2;
    cfg.method   = 'warp';
    elec3 = electrodenormalize(cfg);
    
    transf = elec3.m;
    transfmat = traditional(elec3.m);
    elec1.pnt = transfmat*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
    elec1.pnt = elec1.pnt(1:3,:)';

% redraw GUI
% ----------
function redrawgui(fid)
    dat = get(fid, 'userdata');
    tmpobj = findobj(fid, 'tag', 'pitch'); set(tmpobj, 'string', num2str(dat.transform(4),3));
    tmpobj = findobj(fid, 'tag', 'roll' ); set(tmpobj, 'string', num2str(dat.transform(5),3));
    tmpobj = findobj(fid, 'tag', 'yaw'  ); set(tmpobj, 'string', num2str(dat.transform(6),3));
    tmpobj = findobj(fid, 'tag', 'right'  ); set(tmpobj, 'string', num2str(dat.transform(1),3));
    tmpobj = findobj(fid, 'tag', 'forward'); set(tmpobj, 'string', num2str(dat.transform(2),3));
    tmpobj = findobj(fid, 'tag', 'up'     ); set(tmpobj, 'string', num2str(dat.transform(3),3));
    tmpobj = findobj(fid, 'tag', 'resizex'); set(tmpobj, 'string', num2str(dat.transform(7),3));
    tmpobj = findobj(fid, 'tag', 'resizey'); set(tmpobj, 'string', num2str(dat.transform(8),3));
    tmpobj = findobj(fid, 'tag', 'resizez'); set(tmpobj, 'string', num2str(dat.transform(9),3));
    tmpview = view;
    
    if size(dat.transform,1) > 1
        dat.electransf.pnt = dat.transform*[ dat.elec1.pnt ones(size(dat.elec1.pnt,1),1) ]';
    else
        dat.electransf.pnt = traditional(dat.transform)*[ dat.elec1.pnt ones(size(dat.elec1.pnt,1),1) ]';
    end;
    dat.electransf.pnt = dat.electransf.pnt(1:3,:)';
    set(fid, 'userdata', dat);
    
    h = findobj(fid, 'tag', 'plot3d');
    if isempty(h)
        axis off;
        h = axes('unit', 'normalized', 'position', [0 0.2 1 0.75]);
        set(h, 'tag', 'plot3d');
        axis off;
    else 
        axes(h);
        axis off;
    end;
    plotelec(dat.electransf, 'r');
    if ~isempty(dat.elec2)
        plotelec(dat.elec2     , 'b');
    end;
    set(h, 'tag', 'plot3d');
    
    % plot mesh
    % ---------
    if ~isempty(dat.meshtri) & isempty(findobj(gcf, 'tag', 'mesh'))
        p1 = plotmesh(dat.meshtri, dat.meshpnt, [], 1);
        set(p1, 'tag', 'mesh');
    end;
    
    %view(tmpview);
    rotate3d on    
    
    
    
    
    
    
    
    
  