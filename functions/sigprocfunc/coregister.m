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
%    'chaninfo1' - [struct] channel information structure for first eeglab
%                  channel location structure (might contain fiducials).
%    'chaninfo2' - [struct] channel information structure for second eeglab
%                  channel location structure (might contain fiducials).
%    'autoscale' - ['on'|'off'] autoscale electrode radius when aligning 
%                  fiducials default is 'on'.
%    'helpmsg'   - ['on'|'off'] pop-up help message when calling function.
%                  Default is 'off'.
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
% Revision 1.15  2006/01/13 20:18:04  arno
% invert fiducials etc...
%
% Revision 1.14  2006/01/13 00:45:56  arno
% adding fiducials to the plotted electrodes
%
% Revision 1.13  2006/01/13 00:41:45  arno
% tags etc...
%
% Revision 1.12  2006/01/13 00:34:21  arno
% detect if transformation is failing
%
% Revision 1.11  2006/01/12 23:47:04  arno
% channel conversion to fieldtrip
%
% Revision 1.10  2006/01/12 23:05:08  arno
% allowing fiducials
%
% Revision 1.9  2006/01/11 00:34:36  arno
% noze tag
%
% Revision 1.8  2006/01/11 00:32:55  arno
% edit help message.
%
% Revision 1.7  2006/01/10 23:27:32  arno
% default electrode montages
%
% Revision 1.6  2006/01/10 22:59:11  arno
% move noze
%
% Revision 1.5  2006/01/10 21:54:06  arno
% plotting the noze
%
% Revision 1.4  2006/01/10 00:43:46  arno
% finalizing GUI etc...
%
% Revision 1.3  2005/11/22 00:23:35  arno
% channel correspondance etc...
%
% Revision 1.2  2005/11/08 23:08:55  arno
% header
%
% Revision 1.1  2005/11/08 22:48:05  arno
% Initial revision
%

% Here we bypass channel info
%    'chaninfo'  - [struct] channel information structure for the first
%                  dataset (may use a cell array { struct struct } to enter
%                  channel location info for both channel location struct.
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
        [clist1 clist2] = pop_chancoresp( dat.elec1, dat.elec2, 'autoselect', 'fiducials');
        try,
            [ tmp transform ] = align_fiducials(dat.elec1, dat.elec2, dat.elec2.label(clist2));
            if ~isempty(transform), dat.transform = transform; end;
        catch,
            warndlg2('Transformation failed, try warping fiducials + 1 vertex electrode');
        end;
    elseif strcmpi(com, 'warp')
        [clist1 clist2] = pop_chancoresp( dat.elec1, dat.elec2, 'autoselect', 'all');
        % copy electrode names
        if ~isempty(clist1)
            tmpelec2 = dat.elec2;
            for index = 1:length(clist2)
                tmpelec2.label{clist2(index)} = dat.elec1.label{clist1(index)};
            end;
            [ tmp dat.transform ] = warp_chans(dat.elec1, dat.elec2, tmpelec2.label(clist2));
        end;
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
                            'chaninfo1'  'struct' {}     struct('no', {}); % empty structure
                            'chaninfo2'  'struct' {}     struct('no', {}); % empty structure
                            'transform'  'real'  []      [];
                            'autoscale'  'string' { 'on' 'off' } 'on';
                            'helpmsg'    'string' { 'on' 'off' } 'off';
                            'mesh'       ''      []   defaultmesh });
if isstr(g), error(g); end;

% help message
% ------------
if strcmpi(g.helpmsg, 'on')
    warndlg2( strvcat( 'User channels (sometimes hidden by 3D mesh) are in green, template channels in red.', ...
            'Press ''Warp'' to automatically warp channels to corresponding template channel locations.', ...
            'If desired, then edit the transformation manually.', ...
            ' ', ...
            'To use the location of the corresponding template channels (and discard your current locations),', ...
            'go to menu item "Edit > Channel locations", press the "Look up loc" button, and select a', ...
            'head model. Then re-open DIPFIT "Head model and settings" and select the "No coreg" option.'), ...
          'Co-register channel locations');
end;

% load mesh if any
% ----------------
if ~isempty(g.mesh)
    if isstr(g.mesh)
        try
            g.mesh  = load(g.mesh);
        catch, g.mesh = [];
        end;
    end;
    if ~isempty(g.mesh)
        if isstruct(g.mesh)
            if isfield(g.mesh, 'vol')
                if isfield(g.mesh.vol, 'r')
                    [X Y Z] = sphere(50);
                    dat.meshpnt = { X*max(g.mesh.vol.r) Y*max(g.mesh.vol.r) Z*max(g.mesh.vol.r) };
                    dat.meshtri = [];
                else
                    dat.meshpnt = g.mesh.vol.bnd(1).pnt;
                    dat.meshtri = g.mesh.vol.bnd(1).tri;
                end;
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
else
    dat.meshpnt = [];
    dat.meshtri = [];
end;

% transform to arrays chan1
% -------------------------
TMP                           = eeg_emptyset;
[TMP.chanlocs tmp2 tmp3 ind1] = readlocs(chan1, 'defaultelp', 'besa');
TMP.chaninfo                  = g.chaninfo1;
TMP.nbchan = length(TMP.chanlocs);
cfg   = eeglab2fieldtrip(TMP, 'chanloc_withfid');
elec1 = cfg.elec;

% transform to arrays chan2
% -------------------------
if ~isempty(chan2)
    TMP   = eeg_emptyset;
    [TMP.chanlocs tmp2 tmp3 ind1] = readlocs(chan2, 'defaultelp', 'besa');
    TMP.chaninfo                  = g.chaninfo2;
    TMP.nbchan = length(TMP.chanlocs);
    cfg   = eeglab2fieldtrip(TMP, 'chanloc_withfid');
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

        [ electransf transform ] = align_fiducials(electmp, elec2, g.alignfid);
        if ~isempty(transform), dat.transform = transform; end;
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
dat.elecshow1  = 1:length(elec1.label);
dat.elecshow2  = 1:length(elec2.label);
dat.color1     = [0 1 0];
dat.color2     = [1 0 0];
dat.label1     = 0;
dat.label2     = 0;
dat.meshon     = 1;
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
    h = uicontrol( opt{:}, [0.5  .1  .2  .05], 'style', 'text', 'string', 'Move right {x in mm}');
    h = uicontrol( opt{:}, [0.5  .05 .2  .05], 'style', 'text', 'string', 'Move forward {y in mm}' );
    h = uicontrol( opt{:}, [0.5  0   .2  .05], 'style', 'text', 'string', 'Move up {z in mm}');
    h = uicontrol( opt{:}, [0.7  .1  .1  .05], 'tag', 'right'  , 'callback', cbright  , 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.7  .05 .1  .05], 'tag', 'forward', 'callback', cbforward, 'style', 'edit', 'string', '' );
    h = uicontrol( opt{:}, [0.7  0   .1  .05], 'tag', 'up'     , 'callback', cbup     , 'style', 'edit', 'string', '');
    h = uicontrol( opt{:}, [0.8  .1  .2  .05], 'style', 'pushbutton', 'string', 'Align fiducials', 'callback', cb_fid);
    h = uicontrol( opt{:}, [0.8  .05 .2  .05], 'style', 'pushbutton', 'string', 'Warp montage', 'callback', cb_warp );
    h = uicontrol( opt{:}, [0.8  0   .1  .05], 'style', 'pushbutton', 'string', 'Cancel', 'callback', 'close(gcbf);' );
    h = uicontrol( opt{:}, [0.9  0   .1  .05], 'style', 'pushbutton', 'string', 'Ok', 'tag', 'ok', 'callback', cb_ok);
    
    % put labels next to electrodes
    % -----------------------------
    cb_label1 = [   'tmp = get(gcbf, ''userdata'');' ...
                    'if tmp.label1, set(gcbo, ''string'', ''Labels on'');' ...
                    'else           set(gcbo, ''string'', ''Labels off'');' ...
                    'end;' ...
                    'tmp.label1 = ~tmp.label1;' ...
                    'set(gcbf, ''userdata'', tmp);' ...
                    'clear tmp;' ...
                    'coregister(''redraw'', gcbf);' ];
    cb_label2 = [   'tmp = get(gcbf, ''userdata'');' ...
                    'if tmp.label2, set(gcbo, ''string'', ''Labels on'');' ...
                    'else           set(gcbo, ''string'', ''Labels off'');' ...
                    'end;' ...
                    'tmp.label2 = ~tmp.label2;' ...
                    'set(gcbf, ''userdata'', tmp);' ...
                    'clear tmp;' ...
                    'coregister(''redraw'', gcbf);' ];
    cb_mesh   = [   'tmp = get(gcbf, ''userdata'');' ...
                    'if tmp.meshon, set(gcbo, ''string'', ''Mesh on'');' ...
                    'else           set(gcbo, ''string'', ''Mesh off'');' ...
                    'end;' ...
                    'tmp.meshon = ~tmp.meshon;' ...
                    'set(gcbf, ''userdata'', tmp);' ...
                    'clear tmp;' ...
                    'coregister(''redraw'', gcbf);' ];
    cb_elecshow1 = [ 'tmp = get(gcbf, ''userdata'');' ...
                      'tmp.elecshow1 = pop_chansel( tmp.elec1.label, ''select'', tmp.elecshow1 );' ...
                      'set(gcbf, ''userdata'', tmp);' ...
                      'clear tmp;' ...
                      'coregister(''redraw'', gcbf);' ];
    cb_elecshow2 = [ 'tmp = get(gcbf, ''userdata'');' ...
                      'tmpstrs = { ''19 elec in 10/20'' ''33 elec in 10/20'' ''61 elec in 10/20'' ''all 347 elec'' };' ...
                      'tmpres = inputgui( ''uilist'', {{ ''style'' ''text'' ''string'' ''show only'' } ' ...
                                        ' { ''style'' ''listbox'' ''string'' strvcat(tmpstrs) }}, ' ...
                                        ' ''geometry'', { 1 1 }, ''geomvert'', [1 3] );' ...
                      'if ~isempty(tmpres), tmp.elecshow2 = tmpstrs{tmpres{1}}; end;' ... 
                      'set(gcbf, ''userdata'', tmp);' ...
                      'clear tmp tmpres;' ...
                      'coregister(''redraw'', gcbf);' ];
    h = uicontrol( opt{:}, [0  0.95 .13 .05], 'style', 'pushbutton', 'backgroundcolor', dat.color1, 'string', 'Labels on', 'callback', cb_label1 );
    h = uicontrol( opt{:}, [0  0.9  .13 .05], 'style', 'pushbutton', 'backgroundcolor', dat.color1, 'string', 'Electrodes', 'callback', cb_elecshow1 );

    h = uicontrol( opt{:}, [0 0.85  .13 .05], 'style', 'pushbutton', 'backgroundcolor', dat.color2, 'string', 'Labels on', 'callback', cb_label2 );
    h = uicontrol( opt{:}, [0 0.8   .13 .05], 'style', 'pushbutton', 'backgroundcolor', dat.color2, 'string', 'Electrodes', 'callback', cb_elecshow2 );
    
    h = uicontrol( opt{:}, [0 0.75  .13 .05], 'style', 'pushbutton', 'string', 'Mesh off', 'callback', cb_mesh );
    
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
function plotelec(elec, elecshow, color, tag);
    
    X1 = elec.pnt(elecshow,1);
    Y1 = elec.pnt(elecshow,2);
    Z1 = elec.pnt(elecshow,3);    

    XL = xlim;
    YL = ylim;
    ZL = zlim;
    
    lim=max(1.05*max([X1;Y1;Z1]), max([XL YL ZL]));
    eps=lim/20;
    delete(findobj(gcf, 'tag', tag));
    
    % make bigger if fiducial
    % ------------------------
    fidlist = { 'nz' 'lpa' 'rpa' };
    [tmp fids ] = intersect(lower(elec.label(elecshow)), fidlist);
    nonfids     = setdiff(1:length(elec.label(elecshow)), fids);
    h1 = plot3(X1(nonfids),Y1(nonfids),Z1(nonfids), 'o', 'color', color); hold on;
    set(h1, 'tag', tag, 'marker', '.', 'markersize', 20);
    if ~isempty(fids)
        h2 = plot3(X1(fids),Y1(fids),Z1(fids), 'o', 'color', color*2/3); hold on;
        set(h2, 'tag', tag, 'marker', '.', 'markersize', 35); % make bigger if fiducial
    end;

    % plot axis and labels
    %- -------------------
    if isempty(findobj(gcf, 'tag', 'axlabels'))
        plot3([0.08 0.12],[0 0],[0 0],'r','LineWidth',4) % nose
        plot3([0 lim],[0 0],[0 0],'b--', 'tag', 'axlabels')                 % axis
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
    end;
    axis([-lim lim -lim lim -lim*0.5 lim])
    axis equal;

% decode labels for electrode caps
% --------------------------------
function indices = decodelabels( strchan );
    if ~isstr(strchan), indices = strchan; return; end;
    switch strchan
        case '19 elec in 10/20', indices = [ 1 2 3 4 6 19 21 23 25 27 41 43 45 47 49 63 65 67 69 71 84 86 ];
        case '33 elec in 10/20', indices = [ 1 2 3 4 6 18 19 21 23 25 27 28 31 33 35 37 40 41 43 45 47 49 50 53 55 57 59 62 63 65 67 69 71 72 84 86];
        case '61 elec in 10/20', indices = [ 1 2 3 4 5 6 8 10 12 14 16 19 20 21 22 23 24 25 26 27 31 32 33 34 35 36 37 38 41 42 43 44 45 46 47 48 52 53 54 55 56 57 58 59 60 62 63 64 65 67 68 69 70 71 72 74 76 78 80 82 84 85 86 88];
    otherwise, indices = 1:346;
  end;
    
% plot electrode labels
% ---------------------
function plotlabels(elec, elecshow, color, tag);

    for i = 1:length(elecshow)
        coords = elec.pnt(elecshow(i),:);
        coords = coords*1.07;
        text(coords(1), coords(2), coords(3), elec.label{elecshow(i)},'HorizontalAlignment','center',...
            'VerticalAlignment','middle','Color',color, 'FontSize',10, 'tag', tag)
    end
    
% align fiducials
% ---------------
function [elec1, transf] = align_fiducials(elec1, elec2, fidnames)

    % rename fiducials in template
    % ----------------------------
    fidtemplate = { 'Nz' 'LPA' 'RPA' };
    ind21 = strmatch(fidtemplate{1}, elec2.label, 'exact');
    ind22 = strmatch(fidtemplate{2}, elec2.label, 'exact');
    ind23 = strmatch(fidtemplate{3}, elec2.label, 'exact');
    ind11 = strmatch(fidtemplate{1}, elec1.label, 'exact');
    ind12 = strmatch(fidtemplate{2}, elec1.label, 'exact');
    ind13 = strmatch(fidtemplate{3}, elec1.label, 'exact');
    if isempty(ind11) | isempty(ind12) | isempty(ind13)
        transf = []; return;
    end;
    if isempty(ind21) | isempty(ind22) | isempty(ind23)
        transf = []; return;
    else
         elec2.label{ind21} = fidnames{1};
         elec2.label{ind22} = fidnames{2};
         elec2.label{ind23} = fidnames{3};
    end;
    cfg          = [];
    cfg.elec     = elec1;
    cfg.template = elec2;
    cfg.method   = 'realignfiducial'; 
    cfg.fiducial = fidnames;
    elec3 = electrodenormalize(cfg);
    transf = homogenous2traditional(elec3.m);
    
    % test difference
    % ---------------
    diff1 = mean(mean(abs(elec3.m-traditional(transf))));
    transf(6) = -transf(6);
    diff2 = mean(mean(abs(elec3.m-traditional(transf))));
    if diff1 < diff2, transf(6) = -transf(6); end;
    
    diff1 = mean(mean(abs(elec3.m-traditional(transf))));
    transf(5) = -transf(5);
    diff2 = mean(mean(abs(elec3.m-traditional(transf))));
    if diff1 < diff2, transf(5) = -transf(5); end;
    
    diff1 = mean(mean(abs(elec3.m-traditional(transf))));
    transf(4) = -transf(4);
    diff2 = mean(mean(abs(elec3.m-traditional(transf))));
    if diff1 < diff2, transf(4) = -transf(4); end;

    % rescale if necessary
    % --------------------
    coords2 = elec2.pnt([ind21 ind22 ind23],:); dist_coords2 = sqrt(sum(coords2.^2,2));
    coords1 = elec1.pnt([ind11 ind12 ind13],:); dist_coords1 = sqrt(sum(coords1.^2,2));
    ratio = mean(dist_coords2./dist_coords1);
    transf(7:9) = ratio;
    
    transfmat = traditional(transf);
    elec1.pnt = transfmat*[ elec1.pnt ones(size(elec1.pnt,1),1) ]';
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
    tmpobj = findobj(fid, 'tag', 'pitch'); set(tmpobj, 'string', num2str(dat.transform(4),4));
    tmpobj = findobj(fid, 'tag', 'roll' ); set(tmpobj, 'string', num2str(dat.transform(5),4));
    tmpobj = findobj(fid, 'tag', 'yaw'  ); set(tmpobj, 'string', num2str(dat.transform(6),4));
    tmpobj = findobj(fid, 'tag', 'right'  ); set(tmpobj, 'string', num2str(dat.transform(1),4));
    tmpobj = findobj(fid, 'tag', 'forward'); set(tmpobj, 'string', num2str(dat.transform(2),4));
    tmpobj = findobj(fid, 'tag', 'up'     ); set(tmpobj, 'string', num2str(dat.transform(3),4));
    tmpobj = findobj(fid, 'tag', 'resizex'); set(tmpobj, 'string', num2str(dat.transform(7),4));
    tmpobj = findobj(fid, 'tag', 'resizey'); set(tmpobj, 'string', num2str(dat.transform(8),4));
    tmpobj = findobj(fid, 'tag', 'resizez'); set(tmpobj, 'string', num2str(dat.transform(9),4));
    tmpview = view;
    
    if size(dat.transform,1) > 1
        dat.electransf.pnt = dat.transform*[ dat.elec1.pnt ones(size(dat.elec1.pnt,1),1) ]';
    else
        dat.electransf.pnt = traditional(dat.transform)*[ dat.elec1.pnt ones(size(dat.elec1.pnt,1),1) ]';
    end;
    dat.electransf.pnt   = dat.electransf.pnt(1:3,:)';
    dat.electransf.label = dat.elec1.label;
    set(fid, 'userdata', dat);
    
    h = findobj(fid, 'tag', 'plot3d');
    if isempty(h)
        axis off;
        h = axes('unit', 'normalized', 'position', [0 0.2 1 0.75]);
        set(h, 'tag', 'plot3d');
        axis off;
    else 
        axes(h);
        %axis off;
    end;
    plotelec(dat.electransf, dat.elecshow1, dat.color1, 'elec1');
    if ~isempty(dat.elec2)
        dat.elecshow2 = decodelabels( dat.elecshow2 );
        plotelec(dat.elec2, dat.elecshow2, dat.color2, 'elec2');
    end;
    set(h, 'tag', 'plot3d');
    
    % plot mesh
    % ---------
    if ~isempty(dat.meshpnt) & isempty(findobj(gcf, 'tag', 'mesh'))
        if ~isempty(dat.meshtri)
            p1 = plotmesh(dat.meshtri, dat.meshpnt, [], 1);
            set(p1, 'tag', 'mesh');
        else
            facecolor(1,1,1) = 1; facecolor(1,1,2) = .75; facecolor(1,1,3) = .65;
            cdat = repmat( facecolor, [ size(dat.meshpnt{1}) 1]);
            h = mesh(dat.meshpnt{1}, dat.meshpnt{2}, dat.meshpnt{3}, ...
                'cdata', cdat, 'tag', 'mesh', 'facecolor', squeeze(facecolor), 'edgecolor', 'none');
            hidden off;
            lightangle(45,30);
            lightangle(45+180,30);
            lighting phong
            s = plotnoze([85 0 -75 0 0 pi/2 10 10 40]);
            set(s, 'tag', 'mesh');
        end;
    end;
    meshobj = findobj(gcf, 'tag', 'mesh');
    if dat.meshon
         set( meshobj, 'visible', 'on');
    else set( meshobj, 'visible', 'off');
    end;
    
    % plot electrodes
    % ---------------
    delete(findobj(gcf, 'tag', 'elec1labels'));        
    delete(findobj(gcf, 'tag', 'elec2labels'));        
    if dat.label1
        plotlabels(dat.electransf, dat.elecshow1, dat.color1, 'elec1labels');
    end;
    if dat.label2
        plotlabels(dat.elec2, dat.elecshow2, dat.color2, 'elec2labels');
    end;
    
    %view(tmpview);
    rotate3d on    
    
  
% function to plot the noze
% -------------------------
function s = plotnoze(transf, col)

    if nargin < 1
        transf = [0 0 0 0 0 0 1 1 1];
    end;
    if nargin < 2
        col = [1 0.75 0.65 ];
    end;
        
    x=[ % cube
     NaN -1 1 NaN
      -1 -1 1 1
      -1 -1 1 1
     NaN -1 1 NaN
     NaN -1 1 NaN
     NaN NaN NaN NaN
     ];
    
    y=[ % cube
     NaN -1 -1   NaN
      -1 -1 -1   -1
       1  1  1    1
     NaN  1  1   NaN
     NaN -1 -1   NaN
     NaN NaN NaN NaN
     ];
    
    z=[ % cube
     NaN 0 0 NaN
       0 1 1 0
       0 1 1 0
     NaN 0 0 NaN
     NaN 0 0 NaN
     NaN NaN NaN NaN
     ];

    x=[ % noze
     NaN -1  1 NaN
      -1  0  0 1
     -.3  0  0 .3
     NaN -.3 .3 NaN
     NaN -1  1 NaN
     NaN NaN NaN NaN
     ];
    
    y=[ % noze
     NaN -1 -1   NaN
      -1 -1 -1   -1
       1 -1 -1    1
     NaN  1  1   NaN
     NaN -1 -1   NaN
     NaN NaN NaN NaN
     ];
    
    z=[ % noze
     NaN 0 0 NaN
       0 1 1 0
       0 1 1 0
     NaN 0 0 NaN
     NaN 0 0 NaN
     NaN NaN NaN NaN
     ];

    % apply homogenous transformation
    % -------------------------------
    transfhom = traditional( transf );
    xyz       = [ x(:) y(:) z(:) ones(length(x(:)),1) ];
    xyz2      = transfhom * xyz';
    x(:)      = xyz2(1,:)';
    y(:)      = xyz2(2,:)';
    z(:)      = xyz2(3,:)';

    % dealing with colors
    % -------------------
    cc=zeros(8,3);
    cc(1,:) = col;
    cc(2,:) = col;
    cc(3,:) = col;
    cc(4,:) = col;
    cc(5,:) = col;
    cc(6,:) = col;
    cc(7,:) = col;
    cc(8,:) = col;
    cc(1,:)=[0 0 0]; % black
    cc(2,:)=[1 0 0]; % red
    cc(3,:)=[0 1 0]; % green
    cc(4,:)=[0 0 1]; % blue
    cc(5,:)=[1 0 1]; % magenta
    cc(6,:)=[0 1 1]; % cyan
    cc(7,:)=[1 1 0]; % yellow
    cc(8,:)=[1 1 1]; % white
    cs=size(x);
    c=repmat(zeros(cs),[1 1 3]);
    for i=1:size(cc,1)
        ix=find(x==cc(i,1) &...
                y==cc(i,2) &...
                z==cc(i,3));
        [ir,ic]=ind2sub(cs,ix);
        for k=1:3
            for m=1:length(ir)
                c(ir(m),ic(m),k)=cc(i,k);
            end
        end
    end
    
    % plotting surface
    % ----------------
    facecolor = zeros(size(x,1), size(z,2), 3);
    facecolor(:,:,1) = 1; facecolor(:,:,2) = .75; facecolor(:,:,3) = .65;
   
    s=surf(x,y,z,facecolor);
    set(s, 'edgecolor', [0.5 0.5 0.5]);

    
    
    
    
    
  